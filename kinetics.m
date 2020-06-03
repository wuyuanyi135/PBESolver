function [GD, Bp, Bs] = kinetics(T, sigma, vf)
% GD: Growth (>0) or dissolution (<0) in um/s
% Bp: primary nucleation rate in #/m3/s
% Bs: secondary nucleation rate in #/m3/s
R = 8.3145;

if sigma > 0
%     GD = 1.0e2 * sigma.^(5/6) .* exp(-15000./(R*T));
%     Bp = 8e10 * sigma^(3) .* exp(-35000./(R*T));
%     Bs = 8e15 * sigma^(7/3) .* exp(-35000./(R*T)) .* vf^(0.2);
    GD = 1e-1 * sigma .^ 1;
%     Bp = 0;
%     Bs = sigma.^2 .* (1e8 + 1e12 * (vf/1000) .^ (2/3)); 
    Bp =  sigma.^2 .* 1e8;
    Bs = sigma.^2 .* 1e12 * (vf/1000) .^ (2/3);
else
    Bp = 0;
    Bs = 0;
    % GD = 1.5 * sigma;
    GD = 2.2 * sigma;
end

