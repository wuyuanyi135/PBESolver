function M = particle_moment(lSteps, n, order)
lGrids = (0 : size(n,1)-1)'*lSteps;
M = lSteps .* trapz(n .* lGrids .^ order);
end

