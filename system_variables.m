classdef system_variables
    % Runtime variables that evolve over time.
   
    %% Configuration
    % Cache dependent computation
    properties
        useCache = true
    end
    
    %% State variables governed by the solver
    properties
        concentration
        csd
    end
    
    %% Input variables
    properties
        temperatureC
        inletConcentration
        inletCsd
        residenceTime
    end
    
    %% Dependent properties
    % Do not reuse the system_variable object because the cache will be
    % invalid, unless the useCache has been disabled
    methods
        
        function M = moment(obj, order, dim)
            % Take the first argument as moment order and the second one as
            % integral dimension 
            % TODO: multi-dimension 
            M = memoize(@obj.moment_internal, 'Enabled', obj.useCache);
        end

        function vf = volFrac(obj)
            vf = obj.moment(3, 0) / 1e18 * obj.props.kShape;
        end
    end
    
    methods (Access=private)
        function M = moment_internal(obj, order, dim)
            dL = obj.prop.sizeGrids.interval;
            L = obj.prop.sizeGrids.to_array();
            M = dL * trapz(obj.csd * L^order);
        end
    end
    
    properties(Access=private)
        props system_properties
    end
    
    %% Constructor
    % Take the current property object as argument
    methods
        function obj = system_variables(props)
            obj.props = props;
        end
    end
    
    %% Object creation methods
    methods (Static)
        function obj = create_empty_state(props, initTC)
            % create a system operating at the initial temperature
            obj = system_variables(props);
            obj.temperatureC = initTC;
            obj.concentration = props.solubility(initTC);
            obj.csd = zeros(size(size_grid));
        end
        
        
    end
end

