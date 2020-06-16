function M = particle_moment(lSteps, n, order)
lGrids = (0 : size(n,1)-1)'*lSteps;
Y = n .* (lGrids .^ order);
M = lSteps .* (sum(Y) - 0.5*(Y(1) + Y(end)));
end

