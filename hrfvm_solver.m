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
                f = [0; 0; x; 0]; 
                d = diff(f);
                theta = d(1:end-1)./d(2:end);
                
                % Adopted from CAT.
                theta(isinf(theta)) = 0;
                theta(isnan(theta)) = 2;
                phi = hrfvm_solver.flux_limiter(theta);
                
                if isscalar(GD)
                    % size independent growth
                    x = x - tToL * GD * (f(3:end-1) - f(2:end-2)) ...
                        - tToL * GD / 2 * (1 - tToL * GD) .* ( ...
                        (f(4:end) - f(3:end-1)) .* phi(2:end) - (f(3:end-1) - f(2:end-2)) .* phi(1:end-1) ....
                        );
                else
                    % size dependent growth
                    G = [GD(1); GD(1); GD; GD(end)];
                    Gf = G.*f;
                    x = x - tToL * (Gf(3:end-1) - Gf(2:end-2)) ...
                        - (tToL/2 * G(3:end-1) .* (1 - tToL * G(3:end-1)) .* (f(4:end) - f(3:end-1)) .* phi(1:end-1) ...
                        - tToL/2 * G(2:end-2) .* (1 - tToL * G(2:end-2)) .* (f(3:end-1) - f(2:end-2)) .* phi(2:end));
                end
                % Nucleation after growth, adopted from CAT.
                x(1) = x(1) + (Bs + Bp) * tToL;
            else
                % Dissolving
                % boundary condition
                f = [0; x; 0; 0];
                GD = -GD;
                d = diff(f);
                theta = d(2:end)./d(1:end-1);
                % Adopted from CAT.
                theta(isinf(theta)) = 0;
                theta(isnan(theta)) = 2;
                phi = hrfvm_solver.flux_limiter(theta);
                if isscalar(GD)
                    % size independent dissolution
                    x = x - tToL * GD * (f(2:end-2)- f(3:end-1)) ...
                        - tToL * GD / 2 * (1 - tToL * GD) * ( ...
                        (f(1:end-3)-f(2:end-2)) .* phi(1:end-1) ...
                        - (f(2:end-2)-f(3:end-1)) .* phi(2:end));
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

