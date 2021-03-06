classdef parameterized_system_properties < system_properties
    properties
        %% Solubility definition
        % Polynomial to calculate solubility in kg/kg. The first element is
        % the coefficient of zero-order temperature term.
        solubilityPoly = [4.564e-3; 3.032e-5; 8.437e-6];
        
        %% Kinetics definition
        % Primary nucleation
        pnKp = 1e8
        pnP = 2
        pnKe = 0
        pnEa = 0
        
        % Secondary nucleation
        snKb = 1e10
        snB = 2
        snJ = 2/3
        snEa = 0
        
        % Growth (of each axis)
        % Each element should be of the same size
        gKg = 0.1
        gG = 1
        gBeta = 0
        gEa = 0
        
        % Dissolution (of each axis)
        dKd = 2.2
        dD = 1
        dBeta = 0
        dEa = 0
    end
    methods
        function obj = set.gBeta(obj, x)
            if obj.sizeGrids.highBound * x + 1 <= 0
                error('Invalid gBeta: potential flipped growth rate.')
            end
            obj.gBeta = x;
        end
        function obj = set.dBeta(obj, x)
            if obj.sizeGrids.highBound * x + 1 <= 0
                error('Invalid dBeta: potential flipped dissolution rate.')
            end
            obj.dBeta = x;
        end
    end
    methods
        function obj = parameterized_system_properties()
            obj.sizeGrids = size_grid(999, 999);
        end
    end
    
    methods
        function [GD, Bp, Bs] = kinetics(obj, svar)
            R = 8.3145;
            sigma = svar.sigma;
            vf = svar.vf;
            tC = svar.tC;
            tK = tC + 273.15;
            
            if sigma > 0
                if obj.gBeta == 0
                    % size-independent growth
                    GD = obj.gKg * sigma .^ obj.gG ...
                        * exp(-obj.gEa / R / tK);
                else
                    % size-dependent growth
                    GD = obj.gKg * sigma .^ obj.gG ...
                        * (1 + obj.gBeta * svar.lGrids) ...
                        * exp(-obj.gEa / R / tK);
                end
                Bp =  obj.pnKp * sigma .^ obj.pnP ...
                    * exp(-obj.pnKe / R / log(sigma + 1)^2 ) ...
                    * exp(-obj.pnEa / R / tK);
                
                Bs = obj.snKb * sigma ^ obj.snB * vf ^ obj.snJ * exp(-obj.snEa / R / tK);
            else
                % sigma < 0
                if obj.dBeta == 0
                    % size-independent dissolution
                    GD = obj.dKd * sigma .^ obj.dD ...
                        * exp(-obj.dEa / R / tK);
                else
                    % size-dependent dissolution
                    GD = obj.dKd * sigma .^ obj.dD ...
                        * (1 + obj.dBeta * svar.lGrids) ...
                        * exp(-obj.dEa / R / tK);
                end
                Bp = 0;
                Bs = 0;
            end
        end
        function cStar = solubility(obj, T)
            poly = obj.solubilityPoly;
            Ts = (ones(size(poly)) * T) .^ ((0:numel(poly)-1)');
            cStar = sum(Ts .* poly);
        end
    end
end

