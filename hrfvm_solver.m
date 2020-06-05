classdef hrfvm_solver < pbe_solver
    methods
        function nextState = step(obj, currentState, inputs, cfl, tSpan)
            tNow = 0;
            s = currentState;
            lStep = obj.props.sizeGrids.interval();
            
            cKsDR = obj.props.kShape * obj.props.densityRatio / 1e18;
            x = s.csd;
            c = s.conc;
            m3 = s.moment3;

            if isempty(s.moment3)
                warning('The moment is not set in the initial structure'); 
                m3 = particle_moment(lStep, x, 3);
            end
            
            % Compute solubility
            cStar = obj.props.solubility(inputs.tC);
            while tNow < tSpan
                vf = m3 / 1e18 * obj.props.kShape;
                sigma = c / cStar - 1;
                
                svar = struct();
                svar.tC = inputs.tC;
                svar.vf = vf;
                svar.sigma = sigma;
                % Kinetics
                [GD, Bs, Bp] = obj.props.kinetics(svar);
                
                % early stop
                if max(abs(GD)) <= obj.options.earlyStopThreshold
                    break;
                end

                % estimate time step
                massRate = particle_moment(lStep, x .* GD, 2) * cKsDR;
                massStepLimit = abs((c - cStar) / massRate);
                
                csdTimeLimit = abs(lStep./GD);
                timeLimits = [massStepLimit; csdTimeLimit];
                tStep = min(timeLimits) * cfl;
                tNow = tNow + tStep;
                if tNow > tSpan 
                    tNow = tSpan;
                    tStep = tSpan - tNow;
                end
                
                if sigma > 0
                    % growth
                    % boundary condition
                    x(1) = (Bp + Bs)/GD;
                    f = [x; 0]; 
                    d = diff(f);
                    theta = d(1:end-1)./d(2:end);
                    theta(isnan(theta) | isinf(theta)) = 0;
                    phi = obj.flux_limiter(theta);
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
                    phi = obj.flux_limiter(theta);
                    phinp1 = [phi(1); phi(2:end)];
                    x = x - tStep * GD /lStep * (f(2:end-1)- f(3:end)) ...
                       - (tStep * GD/2/lStep .* (1 - tStep*GD/lStep) .* (f(1:end-2)-f(2:end-1)) .* phi ...
                       - tStep * GD/2/lStep .* (1 - tStep*GD/lStep) .* (f(2:end-1)-f(3:end)) .* phinp1); 
                end
                
                % derive mass change
                newM3 = particle_moment(lStep, x, 3);
                deltaMass = (newM3 - m3) * cKsDR;
                
                m3 = newM3;
                
                c = c - deltaMass;
            end
            
            nextState = struct();
            nextState.csd = x;
            nextState.moment3 = m3;
            nextState.conc = c;
        end
    end
    
    methods (Static)
        function phi = flux_limiter(theta)
            phi = (theta + abs(theta))./(1+abs(theta)); 
        end
    end
end

