classdef system_properties
    % System property such as kinetics, thermodynamics, and etc.
    
    properties
        % Size grid in um. Multi-dimensional size grid is specified in
        % multiple elements of size_grid_spec
        sizeGrids size_grid
        
        % Moment-to-volume shape factor.
        kShape = 0.48
        
        % Density of particle divided by density of solvent
        densityRatio = 1.54
    end
    
    properties
        % Function to calculate solubility kg solid/kg solvent by passing
        % in temperature in degree C
        solubility
        
        % Function to calculate kinetics by passing in the run-time state
        % object
        kinetics
    end
end

