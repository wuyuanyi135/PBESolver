classdef upwind_solver < pbe_solver

methods (Access = protected)
        function x = step_csd(obj, x, Bp, Bs, GD, tStep, lStep)
            tToL = tStep/lStep;
            if GD(1) > 0
                % growth
                % boundary condition
                x(1) = (Bp + Bs)/GD(1);                
                
                if isscalar(GD)
                    % size independent growth
                    x(2:end) = x(2:end) - tToL * GD * (x(2:end) - x(1:end-1));
                else
                    % size dependent growth
                    Gx = GD.*x;
                    x(2:end) = x(2:end) - tToL * (Gx(2:end) - Gx(1:end-1));
                end
            else
                % Dissolving
                % boundary condition
                f = [x; 0];
                GD = -GD;
                if isscalar(GD)
                    % size independent dissolution
                    x = x - tToL * GD  * (f(1:end-1)- f(2:end));
                else
                    % size dependent dissolution
                    D = [GD; 0];
                    Df = D .* f;
                    x = x - tToL * (Df(1:end-1)- Df(2:end));
                end
            end 
        end
    end
end

