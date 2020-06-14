function M = particle_moment(lSteps, n, order)
lGrids = (0 : lSteps : size(n,1))';
M = lSteps .* trapz(n .* lGrids .^ order);
end

