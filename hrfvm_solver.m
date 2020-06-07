classdef hrfvm_solver < pbe_solver %#codegen
    methods
        function obj = hrfvm_solver(options, props)
            obj = obj@pbe_solver(options, props);
        end
    end
    methods
        function nextState = step(obj, currentState, inputs, cfl, tSpan)
            % TODO size dependent GD is not allowed!
            if isempty(currentState)
                st = obj.prevStates;
                stateful = true;
            else
                st = currentState;
                stateful = false;
            end
            
            props = obj.props;
            nProps = numel(props);
            assert(numel(st) == nProps);
            
            tNow = 0;
            s = st;
            sizeGrids = [obj.props.sizeGrids];
            lStep = sizeGrids.interval();
            lGrids = cell(1, nProps);
            for i = 1 : nProps
                lGrids{i} = sizeGrids(i).to_array();
            end
            
            kShapes = [props.kShape];
            cKsDR = kShapes .* [props.densityRatio] / 1e18;
            x = cell(size(props));
            for i = 1 : nProps
                x{i} = s(i).csd;
            end
            c = s(1).conc; % only first concentration is used.
            m3 = [s.moment3];
            
            % Compute solubility
            cStars = zeros(size(props));
            for i = 1 : nProps
                cStars(i) = props(i).solubility(inputs.tC);
            end
            
            % Loop variables initialization
            svar = struct();
            GD = zeros(size(props));
            Bs = zeros(size(props));
            Bp = zeros(size(props));
            massRate = zeros(size(props));
            newM3 = zeros(size(props));
            
            while tNow < tSpan
                vf = m3 / 1e18 .* kShapes;
                sigma = c ./ cStars - 1;
                
                for i = 1 : nProps
                    svar.tC = inputs.tC;
                    svar.vf = vf(i);
                    svar.sigma = sigma(i);
                    % Kinetics
                    [GD(i), Bs(i), Bp(i)] = props(i).kinetics(svar);
                end
                
                % early stop;
                if max(abs(GD)) <= obj.options.earlyStopThreshold
                    break;
                end

                % estimate time step
                for i = 1 : nProps
                    massRate(i) = particle_moment(lStep(i), lGrids{i}, x{i} .* GD(i), 2) .* cKsDR(i);
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
                    x{i} = obj.forward_csd(x{i}, sigma(i), Bp(i), Bs(i), GD(i), tStep, lStep(i));
                    % Derive mass change
                    newM3(i) = particle_moment(lStep(i), lGrids{i}, x{i}, 3);
                end
                deltaMass = (newM3 - m3) .* cKsDR;
                
                m3 = newM3;
                
                c = c - sum(deltaMass);
            end
            
            % TODO this order affects code generation
            nextState = repmat(struct('conc', 0, 'moment3', 0, 'csd', zeros(size(x{1}))), 1, nProps);
            for i = 1 : nProps
                nextState(i).csd = x{i};
                nextState(i).moment3 = m3(i);
                nextState(i).conc = c;
            end
            
            if stateful
                obj.prevStates = nextState;
            end
        end
    end
    
    methods (Static)
        function phi = flux_limiter(theta)
            phi = (theta + abs(theta))./(1+abs(theta)); 
        end
        
        function x = forward_csd(x, sigma, Bp, Bs, GD, tStep, lStep)
            if sigma > 0
                % growth
                % boundary condition
                x(1) = (Bp + Bs)/GD;
                f = [x; 0]; 
                d = diff(f);
                theta = d(1:end-1)./d(2:end);
                theta(isnan(theta) | isinf(theta)) = 0;
                phi = hrfvm_solver.flux_limiter(theta);
                phin_1 = [phi(2:end); phi(end)];

                x(2:end) = x(2:end) - tStep * GD /lStep * (f(2:end-1) - f(1:end-2)) ...
                   - (tStep * GD/2/lStep .* (1 - tStep*GD/lStep) .* (f(3:end) - f(2:end-1)) .* phi ...
                   - tStep * GD/2/lStep .* (1 - tStep*GD/lStep) .* (f(2:end-1) - f(1:end-2)) .* phin_1); 
            else
                % Dissolving
                % boundary condition
                f = [0; x; 0];
                d = diff(f);
                GD = -GD;
                theta = d(2:end)./d(1:end-1);
                theta(isnan(theta) | isinf(theta)) = 0;
                phi = hrfvm_solver.flux_limiter(theta);
                phinp1 = [phi(1); phi(2:end)];
                x = x - tStep * GD /lStep * (f(2:end-1)- f(3:end)) ...
                   - (tStep * GD/2/lStep .* (1 - tStep*GD/lStep) .* (f(1:end-2)-f(2:end-1)) .* phi ...
                   - tStep * GD/2/lStep .* (1 - tStep*GD/lStep) .* (f(2:end-1)-f(3:end)) .* phinp1); 
            end 
        end
    end
end

