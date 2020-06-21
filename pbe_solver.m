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
        
        function nextState = step(obj, currentState, inputs, tSpan)
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
            szProps = size(properties);
            assert(numel(s) == nProps);
            
            tNow = 0;
            timeStepScale = obj.options.timeStepScale;
            isMSMPR = obj.options.isMSMPR;
            if isMSMPR
                % Populate continuous inputs
                resTime = inputs.resTime;
                inConc = inputs.inConc;
                inCSDs = inputs.inCSDs;
                residenceTimeStepScale = obj.options.residenceTimeStepScale;
            end
            
            sizeGrids = [properties.sizeGrids];
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
            cStars = zeros(szProps);
            for i = 1 : nProps
                cStars(i) = properties(i).solubility(inputs.tC);
            end
            
            % Loop variables initialization
            svar = struct();
            GD = cell(szProps);
            Bs = zeros(szProps);
            Bp = zeros(szProps);
            massRate = zeros(szProps);
            newM3 = zeros(szProps);
            
            %% Locate the current in-use size range and Step-artifact prevention
            % This optimization should help reducing the computational load
            % by excluding the unused channels in integration.
            % When continuous crystallization is involved, the effective
            % CSD is the union of the state CSD and the inlet CSD.
            optUseSubCSD = obj.options.useSubCSD;
            effectiveCSDs = x;
            effectiveLGrids = lGrids;
            if optUseSubCSD

                rightSizeCaps = zeros(szProps);
                rPositions = zeros(szProps);
            
                for i = 1 : nProps
                    mark = find(x{i} > eps , 1, 'last');
                    if isempty(mark)
                        mark = 1; % first channel must be zero (nucleation point)
                    end
                    if isMSMPR
                        inCSDMark = find(inCSDs{i} > eps , 1, 'last');
                        if isempty(inCSDMark)
                            inCSDMark = 1;
                        end
                        mark = max(mark, inCSDMark);
                    end
                    rightSizeCaps(i) = lGrids{i}(mark);
                    % Find the first channel that is above the current cap.
                    % rPositions(i) = find(lGrids{i}>=cap, 1);
                    rPositions(i) = mark;
                    % Because each step the growth/dissolution cannot
                    % propagate more than one channel, leaving one spare is
                    % enough.
                    if rPositions(i) == numel(lGrids{i})
                        idx = 1:rPositions(i);
                    else
                        idx = 1:rPositions(i)+1;
                    end
                    effectiveLGrids{i} = lGrids{i}(idx);
                    effectiveCSDs{i} = x{i}(idx);
                end
            end
            
            %% Loop
            while tNow < tSpan
                vf = m3 / 1e18 .* kShapes;
                sigma = c ./ cStars - 1;
                
                for i = 1 : nProps
                    svar.tC = inputs.tC;
                    svar.vf = vf(i);
                    svar.sigma = sigma(i);
                    svar.lGrids = effectiveLGrids{i};
                    % Kinetics
                    [GD{i}, Bs(i), Bp(i)] = properties(i).kinetics(svar);
                end
                
                % early stop
                if obj.check_early_stop(tNow, tSpan, GD, sigma)
                    break;
                end
                
                % calculate time step with the predicted change of
                %   1. concentration change (should not bypass the
                %   solubility line)
                %   2. CSD change (should respect stability condition)
                % TODO: get_time_step slot
                % TODO: when no crystal and dissolving, remove the time
                % step limit on PBE.
                for i = 1 : nProps
                    massRate(i) = particle_moment(lStep(i), effectiveCSDs{i} .* GD{i}, 2) .* cKsDR(i);
                end
                
                massStepLimit = abs((c - cStars) ./ massRate);
                
                csdTimeLimits = zeros(size(GD));
                % TODO size_dependent_mismatch: when disabling SubCSD optimization, and if GD is
                % increasing over size, the time step of unused channels
                % will result in smaller time steps and affects the result
                % consistency between SubCSD ON and SubCSD OFF.
                for i = 1 : numel(GD)
                    if optUseSubCSD
                        noParticle = rPositions(i) == 1;
                    else
                        noParticle = all(effectiveCSDs{i}==0);
                    end
                    if all(GD{i} == 0) || (noParticle && GD{i}(1) < 0)
                        % GD is not useful when no growth/dissolution
                        % or when there is no crystal and it is
                        % dissolving.
                        csdTimeLimits(i) = inf;
                    else
                        csdTimeLimits(i) = min(abs(lStep ./ GD{i}));
                    end   
                end
                csdTimeLimit = min(csdTimeLimits);
                
                % TODO: for MSMPR, the time step and tau should be related
                if isMSMPR
                    resTimeStep = resTime * residenceTimeStepScale;
                    timeLimits = [massStepLimit csdTimeLimit resTimeStep];
                else
                    timeLimits = [massStepLimit csdTimeLimit];
                end
                
                tStep = min(timeLimits) * timeStepScale;
                
                tNow = tNow + tStep;
                if tNow > tSpan 
                    tNow = tSpan;
                    tStep = tSpan - tNow;
                end
                
                for i = 1 : nProps
                    if isMSMPR
                        inCSD = inCSDs{i};
                        n = effectiveCSDs{i};

                        if inCSD == 0
                            % Inlet is clear.
                            effectiveCSDs{i} = n + tStep * (0 - n)/resTime;
                        else
                            % Inlet contains crystals.
                            if optUseSubCSD
                                % inCSD may have different effective CSD range.
                                % In this branch, the current CSD and the
                                % effective CSD range should be merged.
                                nIn = inCSD(1: numel(n));
                            else
                                % Assume that the inlet CSD should always be
                                % the same size as the effective CSD.
                                nIn = inCSD;
                            end
                            effectiveCSDs{i} = n + tStep * (nIn - n)/resTime;
                        end 
                    end
                    
                    if optUseSubCSD
                        gd = GD{i};
                        % Determine next CSD right end location
                        rSzCap = rightSizeCaps(i);
                        if isscalar(gd)
                            rightSizeCaps(i) = rSzCap + gd * tStep;
                        else
                            gdInterp = interp1(effectiveLGrids{i}, gd, rSzCap);
                            rightSizeCaps(i) = rSzCap + gdInterp * tStep;
                        end

                        % After update, the right size cap may enter next
                        % channel. That's why one extra channel is included in 
                        % the above code. If the size cap does not move to next
                        % channel, the extra calculated channel should be
                        % discarded.
                        oldPosition = rPositions(i);
                        newPosition = find(lGrids{i}>=rightSizeCaps(i), 1);
                        if isempty(newPosition)
                            % All channels are overflown
                            newPosition = numel(lGrids{i});
                        end
                        effCSD = effectiveCSDs{i};
                        % Calculate next CSD
                        effCSD = obj.step_csd(effCSD, Bp(i), Bs(i), gd, tStep, lStep(i));

                        if newPosition == oldPosition
                            % Growing, dissolving, anything, but not
                            % entering other channel. oldPosition is the
                            % previous right side position. We use it to
                            % drop the right extra channel.
                            if numel(effCSD) < numel(lGrids{i})
                                effCSD(end) = 0;
                            end
                            effectiveCSDs{i} = effCSD;

                        elseif newPosition == oldPosition - 1
                            % Dissolution-induced mismatch. Drop the old
                            % channel
                            effectiveCSDs{i} = effCSD(1:newPosition);
                            % Update the right (old) position & effective LGrids
                            rPositions(i) = newPosition;
                            idx = 1:(newPosition+1);
                            effectiveLGrids{i} = lGrids{i}(idx);
                        elseif newPosition == oldPosition + 1
                            % Growth-induced mismatch, use full channels
                            effectiveCSDs{i} = [effCSD; 0];
                            % Update the right (old) position & effective LGrids
                            rPositions(i) = newPosition;
                            
                            if newPosition == numel(lGrids{i})
                                idx = 1:(newPosition);
                            else
                                idx = 1:(newPosition+1);
                            end
                            
                            effectiveLGrids{i} = lGrids{i}(idx);
                        else
                            error('Unexpected branch reached.')
                        end
                    else
                        effectiveCSDs{i} = obj.step_csd(effectiveCSDs{i}, Bp(i), Bs(i), GD{i}, tStep, lStep(i));
                    end
                    
                    % Derive mass change
                    newM3(i) = particle_moment(lStep(i), effectiveCSDs{i}, 3);
                end
                deltaMass = (newM3 - m3) .* cKsDR;
                
                m3 = newM3;
                
                c = c - sum(deltaMass);
                if isMSMPR
                    % TODO use higher order solver scheme.
                    c = c + tStep * (inConc - c) / resTime;
                end
            end
            
            % TODO this order affects code generation
            nextState = repmat(struct('conc', 0, 'moment3', 0, 'csd', zeros(size(x{1}))), 1, nProps);
            for i = 1 : nProps
                effCSD = effectiveCSDs{i};
                csd = zeros(size(lGrids{i}));
                csd(1:numel(effCSD)) = effCSD;
                nextState(i).csd = csd;
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
            GDs = zeros(size(GD));
            for i = 1 : numel(GD)
                GDs(i) = max(abs(GD{i}));
            end
            
            shouldEarlyStop = max(GDs) < obj.options.earlyStopThreshold;
        end
      
    end
end

