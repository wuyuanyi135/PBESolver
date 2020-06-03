function [xs, cs, ts, GDs, Bss, sigmas] = hrfvm_solver( x0, c0, cfl, tspan, T, lstep, k_s, density_ratio, crossing_early_stop_speed) %#codegen
% A fixed time step, 1d high resolution fvm implementation for PBE solution
% Breakage and agglomeration is not supported.
% If sigma is zero, the solver will quit.
% This solver is used for batch system. Assume temperature is a input that
% is well controlled by an external system.
% x0: 1xN current number density
% T: temperature in kelvin
% density_ratio: rho_p / rho_solvent
% k_s: shape factor k_s * L^3 = actual volume
% x: return value t_point x num_size_grid. The first point is not included.
% c: return value of concentration t_point-sized vector
% ts: relative time points

xs = [];
cs = [];
ts = [];

x = x0;
c = c0;
t = 0;
L = (0:(size(x0)-1))' * lstep;
% In one calculation loop temperature does not change.
c_star = solubility(T);

GDs = [];
Bss = [];
sigmas = [];
while t < tspan
    % get solubility
    sigma = c/c_star - 1;
    
    % calculate current volume fraction
    M3_current =  particle_moment(lstep, x, 3)/ 1e18;
    v_f = M3_current * k_s ;
    
    % compute kinetics
    [GD, Bp, Bs] = kinetics(T, sigma, v_f);
    % TODO
%     GDs = [GDs; GD];
%     Bss = [Bss; Bs];
%     sigmas = [sigmas; sigma];
    % early stop
    if abs(GD) <= crossing_early_stop_speed % todo
%         ts = [ts; tspan];
%         cs = [cs; cs(end)];
%         xs = [xs; xs(end, :)];
        return;
    end
    
    mass_rate = particle_moment(lstep, x, 2) * GD / 1e18 * k_s * density_ratio;
    conc_tstep_limit = abs((c - c_star) / mass_rate);
    tstep = min([lstep./abs(GD)*cfl; conc_tstep_limit*cfl]);
    
    t_new = t + tstep;
    if t_new > tspan
        t_new = tspan;
        tstep = tspan - t;
    end
    
    % derive next step of csd with hr fvm
    if sigma > 0
        % growth
        % boundary condition
        x(1) = (Bp + Bs)/GD;
        f = [x; 0]; 
        d = diff(f);
        theta = d(1:end-1)./d(2:end);
        theta(isnan(theta) | isinf(theta)) = 0;
        phi = flux_limiter(theta);
        phin_1 = [phi(2:end); phi(end)];
        x_new = x;
        x_new(2:end) = x(2:end) - tstep * GD /lstep * (f(2:end-1) - f(1:end-2)) ...
           - (tstep * GD/2/lstep .* (1 - tstep*GD/lstep) .* (f(3:end) - f(2:end-1)) .* phi ...
           - tstep * GD/2/lstep .* (1 - tstep*GD/lstep) .* (f(2:end-1) - f(1:end-2)) .* phin_1); 
       %x_new(1) = f(1); %TODO
    else
        % Dissolving
        % boundary condition
        f = [0; x; 0];
        d = diff(f);
        GD = -GD;
        theta = d(2:end)./d(1:end-1);
        theta(isnan(theta) | isinf(theta)) = 0;
        phi = flux_limiter(theta);
        phinp1 = [phi(1); phi(2:end)];
        x_new = x - tstep * GD /lstep * (f(2: end-1)- f(3:end)) ...
           - (tstep * GD/2/lstep .* (1 - tstep*GD/lstep) .* (f(1:end-2)-f(2:end-1)) .* phi ...
           - tstep * GD/2/lstep .* (1 - tstep*GD/lstep) .* (f(2:end-1)-f(3:end)) .* phinp1); 
    end
    
    % calculate delta mass (kg) / kg solution
    delta_mass = (particle_moment(lstep, x_new, 3)/1e18 - M3_current) * k_s * density_ratio;
    
    c_new = c - delta_mass;
    
%     ts = [ts; t_new];
%     cs = [cs; c_new];
%     xs = [xs; x_new'];
    
    x = x_new;
    c = c_new;
    t = t_new;
end

end

