classdef hrfvm_solver < pbe_solver %#codegen
    methods
        function obj = hrfvm_solver(options, props)
            obj = obj@pbe_solver(options, props);
        end
    end
    
    methods (Static)
        function phi = flux_limiter(theta)
            phi = (theta + abs(theta))./(1+abs(theta)); 
        end
    end
    methods (Access = protected)
        function x = step_csd(obj, x, Bp, Bs, GD, tStep, lStep)
            tToL = tStep/lStep;
            if GD(1) > 0
                % growth
                % boundary condition
                x(1) = (Bp + Bs)/GD(1);
                f = [x; 0]; 
                d = diff(f);
                theta = d(1:end-1)./d(2:end);
                theta(isnan(theta) | isinf(theta)) = 0;
                phi = hrfvm_solver.flux_limiter(theta);
                phin_1 = [phi(2:end); phi(end)];
                
                % Note, during growth, the first channel is governed by the
                % boundary condition. Hence, we do not prepend annother
                % channel to f. However, for f_n or G_n in the equation,
                % they should be f(2:end-1) or G(2:end-1).
                
                if isscalar(GD)
                    % size independent growth
                    x(2:end) = x(2:end) - tStep * GD /lStep * (f(2:end-1) - f(1:end-2)) ...
                        - (tStep * GD/2/lStep .* (1 - tStep*GD/lStep) .* (f(3:end) - f(2:end-1)) .* phi ...
                        - tStep * GD/2/lStep .* (1 - tStep*GD/lStep) .* (f(2:end-1) - f(1:end-2)) .* phin_1);
                else
                    % size dependent growth
                    G = [GD; GD(end)];
                    Gf = G.*f;
                    x(2:end) = x(2:end) - tToL * (Gf(2:end-1) - Gf(1:end-2)) ...
                        - (tToL/2 * G(2:end-1) .* (1 - tToL * G(2:end-1)) .* (f(3:end) - f(2:end-1)) .* phi ...
                        - tToL/2 * G(1:end-2) .* (1 - tToL * G(2:end-1)) .* (f(2:end-1) - f(1:end-2)) .* phin_1);
                end
            else
                % Dissolving
                % boundary condition
                f = [0; x; 0];
                GD = -GD;
                d = diff(f);
                theta = d(2:end)./d(1:end-1);
                theta(isnan(theta) | isinf(theta)) = 0;
                phi = hrfvm_solver.flux_limiter(theta);
                phinp1 = [phi(1); phi(2:end)];
                if isscalar(GD)
                    % size independent dissolution
                    x = x - tStep * GD /lStep * (f(2:end-1)- f(3:end)) ...
                        - (tStep * GD/2/lStep .* (1 - tStep*GD/lStep) .* (f(1:end-2)-f(2:end-1)) .* phi ...
                        - tStep * GD/2/lStep .* (1 - tStep*GD/lStep) .* (f(2:end-1)-f(3:end)) .* phinp1);
                else
                    % size dependent dissolution
                    D = [GD(1); GD; GD(end)];
                    Gf = D .* f;
                    x = x - tToL * (Gf(2:end-1)- Gf(3:end)) ...
                        - (tToL * D(2:end-1)/2 .* (1 - tToL * D(2:end-1)) .* (f(1:end-2)-f(2:end-1)) .* phi ...
                        - tToL * D(3:end)/2 .* (1 - tToL*D(3:end)) .* (f(2:end-1)-f(3:end)) .* phinp1);
                end
            end 
        end
    end
end

