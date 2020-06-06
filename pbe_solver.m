classdef pbe_solver < handle
    % Base solver
    
    properties
        props
        options
    end
    
    properties
        prevStates
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
        function obj = set_initial_states(obj, states)
            obj.prevStates = states;
        end
        
        function nextState = step(currentState, inputs, cfl, tSpan)
            % currentState: structure:
            %   conc: double scalar
            %   csd: double matrix
            %   moment3: double (optional)
            % when the currentState is passed as an array, only the first
            % concentration is respected as the initial concentration. The
            % element number should match the length of the props. They are
            % treated as different polymorphism so that the concentration
            % should be shared.
            % 
            % When the currentState is [], the stored `prevStates` is used
            % and updated for the purpose of stateful solver.
            %
            % inputs: structure:
            %   tC: double - temperature in degC
            %   inletCsd: double - inlet CSD
            %   residenceTime: double
            %   inletConc: double
        end
    end
end

