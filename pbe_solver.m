classdef pbe_solver < handle
    % Base solver
    
    properties
        props
        options
    end
    
    properties
        prevStates
    end
    
    methods
        function obj = pbe_solver(options, props)
            % When props has more than one elements, it represents a
            % polymorphism system.

            obj.options = options;
            obj.props = props;
        end
    end
    
    methods
        function obj = set_initial_states(obj, states)
            obj.prevStates = states;
        end
        
        function nextState = step(obj, currentState, inputs, cfl, tSpan)
            %% Derive next state
            % currentState: structure:
            %   conc: double scalar
            %   csd: double matrix
            %   moment3: double (optional)
            % when the currentState is passed as an array, only the first
            % concentration is respected as the initial concentration. The
            % element number should match the length of the props. They are
            % treated as different polymorphism so that the concentration
            % should be shared.
            % 
            % When the currentState is [], the stored `prevStates` is used
            % and updated for the purpose of stateful solver.
            %
            % inputs: structure:
            %   tC: double - temperature in degC
            %   inletCsd: double - inlet CSD
            %   residenceTime: double
            %   inletConc: double
            
            %% Initialization
            % determine whether the solver is currently stateful or
            % stateless.
            if isempty(currentState)
                s = obj.prevStates;
                stateful = true;
            else
                s = currentState;
                stateful = false;
            end
            
            properties = obj.props;
            nProps = numel(properties);
            assert(numel(s) == nProps);
            
            tNow = 0;
            sizeGrids = [obj.props.sizeGrids];
            lStep = sizeGrids.interval();
            lGrids = cell(1, nProps);
            for i = 1 : nProps
                lGrids{i} = sizeGrids(i).to_array();
            end
            
            kShapes = [properties.kShape];
            cKsDR = kShapes .* [properties.densityRatio] / 1e18;
            x = cell(size(properties));
            for i = 1 : nProps
                x{i} = s(i).csd;
            end
            c = s(1).conc; % only first concentration is used.
            m3 = [s.moment3];
            
            % Compute solubility
            cStars = zeros(size(properties));
            for i = 1 : nProps
                cStars(i) = properties(i).solubility(inputs.tC);
            end
            
            % Loop variables initialization
            svar = struct();
            GD = zeros(size(properties));
            Bs = zeros(size(properties));
            Bp = zeros(size(properties));
            massRate = zeros(size(properties));
            newM3 = zeros(size(properties));
            
            %% Locate the current in-use size range
            % This optimization should help reducing the computational load
            % by excluding the unused channels in integration.
            effectiveCSDs = x;
            effectiveLGrids = lGrids;
            optUseSubCSD = obj.options.useSubCSD;
            if optUseSubCSD
                csdNonZeroMarks = zeros(size(x));
                for i = 1 : nProps
                    mark = find(x{i} ~= 0, 1, 'last');
                    if isempty(mark)
                        mark = 1;
                    end
                    csdNonZeroMarks(i) = mark;
                end
            end
            
            
            %% Loop through
            while tNow < tSpan
                vf = m3 / 1e18 .* kShapes;
                sigma = c ./ cStars - 1;
                
                for i = 1 : nProps
                    svar.tC = inputs.tC;
                    svar.vf = vf(i);
                    svar.sigma = sigma(i);
                    % Kinetics
                    [GD(i), Bs(i), Bp(i)] = properties(i).kinetics(svar);
                end
                
                % early stop
                if obj.check_early_stop(tNow, tSpan, GD, sigma)
                    break;
                end

                % Populate effective CSD (nonzeros + 1 channel)
                if optUseSubCSD
                    for i = 1 : nProps
                        mark = csdNonZeroMarks(i);
                        effectiveCSDs{i} = x{i}(1:mark+1);
                        effectiveLGrids{i} = lGrids{i}(1:mark+1);
                    end
                end
                
                % calculate time step with the predicted change of
                %   1. concentration change (should not bypass the
                %   solubility line)
                %   2. CSD change (should respect stability condition)
                % TODO: get_time_step slot
                for i = 1 : nProps
                    massRate(i) = particle_moment(lStep(i), effectiveLGrids{i}, effectiveCSDs{i} .* GD(i), 2) .* cKsDR(i);
                end
                
                massStepLimit = abs((c - cStars) ./ massRate);
                
                csdTimeLimit = abs(lStep ./ GD);
                timeLimits = [massStepLimit csdTimeLimit];
                tStep = min(timeLimits) * cfl;
                
                tNow = tNow + tStep;
                if tNow > tSpan 
                    tNow = tSpan;
                    tStep = tSpan - tNow;
                end
                
                for i = 1 : nProps
                    % Calculate next CSD
                    effectiveCSDs{i} = obj.step_csd(effectiveCSDs{i}, Bp(i), Bs(i), GD(i), tStep, lStep(i));
                    % Derive mass change
                    newM3(i) = particle_moment(lStep(i), effectiveLGrids{i}, effectiveCSDs{i}, 3);
                end
                deltaMass = (newM3 - m3) .* cKsDR;
                
                m3 = newM3;
                
                c = c - sum(deltaMass);
                
                % Update effective CSD mark
                if optUseSubCSD
                    for i = 1 : nProps
                        % Overwrite the original CSD
                        x{i}(1:csdNonZeroMarks(i)+1) = effectiveCSDs{i};
                        if effectiveCSDs{i}(end) == 0
                            % The extra channel does not have data. Either not
                            % growing or dissolving. Then, the mark should
                            % backtrace the new zero channels
                            mark = find(effectiveCSDs{i} ~= 0, 1, 'last');
                            if isempty(mark)
                                mark = 1;
                            end
                            csdNonZeroMarks(i) = mark;
                        else
                            % Otherwise, right shift 1 channel for next growth.
                            if csdNonZeroMarks(i) < size(x{i}, 1)
                                csdNonZeroMarks(i) = csdNonZeroMarks(i) + 1;
                            end
                        end
                    end
                end
            end
            
            % TODO this order affects code generation
            nextState = repmat(struct('conc', 0, 'moment3', 0, 'csd', zeros(size(x{1}))), 1, nProps);
            for i = 1 : nProps
                if optUseSubCSD
                    nextState(i).csd = x{i};
                else
                    % when not using sub-CSD, the x does not update in the
                    % loop.
                    nextState(i).csd = effectiveCSDs{i};
                end
                nextState(i).moment3 = m3(i);
                nextState(i).conc = c;
            end
            
            if stateful
                obj.prevStates = nextState;
            end
        end
    end
    
    methods (Access = protected)
        function nextCSD = step_csd(obj, x, Bp, Bs, GD, tStep, lStep)
            %% Derive next CSD based on current conditions
        end
        function shouldEarlyStop = check_early_stop(obj, tNow, tSpan, GD, sigma)
            %% Check whether current condition should trigger early stop
            % tNow: cuurrent time
            % tSpan: current solution time span
            % GD: growth or dissolution rates
            % sigma: supersaturations
            shouldEarlyStop = max(abs(GD)) <= obj.options.earlyStopThreshold;
        end
      
    end
end

