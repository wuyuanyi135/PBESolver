function M = particle_moment(lSteps, n, order)
% TODO: If the truncated CSD (n) is passed. Beware that the ending zero
% must be present otherwise the result will be inaccurate.
lGrids = (0 : size(n,1)-1)'*lSteps;
Y = n .* (lGrids .^ order);
% M = lSteps .* (sum(Y) - 0.5*(Y(1) + Y(end)));
M = lSteps .* (sum(Y) - 0.5*(Y(1)));
end

