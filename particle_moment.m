function M = particle_moment(lSteps, lGrids, n, order)
M = lSteps .* trapz(n .* lGrids .^ order);
end

