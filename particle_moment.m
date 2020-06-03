function M = particle_moment(lstep, n, order) %#codegen
L = (lstep * (0: numel(n)-1))'; 
M = trapz(L, n .* L .^ order);
end

