function phi = flux_limiter(theta)
phi = (theta + abs(theta))./(1+abs(theta)); 
end

