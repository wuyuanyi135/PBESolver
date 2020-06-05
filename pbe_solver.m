classdef pbe_solver
    % Base solver
    
    properties
        props
        options
    end
    
    methods
        function obj = pbe_solver(options, props)
            % When props has more than one elements, it represents a
            % polymorphism system.
            obj.options = options;
            obj.props = props;
        end
    end
    
    methods
        function nextState = step(currentState, inputs, cfl, tSpan)
            % currentState: structure:
            %   conc: double scalar
            %   csd: double matrix
            %   moment3: double (optional)
            %
            % inputs: structure:
            %   tC: double - temperature in degC
            %   inletCsd: double - inlet CSD
            %   residenceTime: double
            %   inletConc: double
        end
    end
end

