function c = solubility(T)
% T: temperature in kelvin
% return mass solubility (kg solute/kg solvent)
T = T - 273.15;
%c = (7.644e-3.*T.^2 -1.165e-1 .* T + 6.222)/1000;
c = (8.437e-3 * T.^2 + 3.032e-2 * T + 4.564) / 1000;
end

