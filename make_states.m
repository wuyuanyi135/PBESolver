function states = make_states(conc, initCsd, props)
states = struct();
states.conc = conc;
grid = props.sizeGrids.to_array();
states.moment3 = 0;
if initCsd == 0
    coder.varsize('csd');
    csd = zeros(size(grid));
    states.csd = csd;
    return;
end

if size(grid) == size(initCsd)
    states.csd = initCsd; 
    lStep = props.sizeGrids.interval();
    states.moment3 = particle_moment(lStep, initCsd, 3);
else
    error('initCsd should be 0 or a vector of the same size of sizeGrids')
end
end

