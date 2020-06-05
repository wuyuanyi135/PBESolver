function states = make_states(conc, initCsd, props)
states = struct();
states.conc = conc;
grid = props.sizeGrids.to_array();
states.moment3 = 0;
if initCsd == 0
    csd = zeros(size(grid));
    states.csd = csd;
    return;
end

if size(grid) == size(initCsd)
    states.csd = initCsd; 
    lstep = props.sizeGrids.interval();
    states.moment3 = particle_moment(lstep, csd, 3);
else
    error('initCsd should be a scalar or a vector of the same size of sizeGrids')
end
end

