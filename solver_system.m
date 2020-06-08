classdef solver_system < matlab.System & matlab.system.mixin.Propagates
    properties(Nontunable)
        sampleTimeSecond = 30;
        cfl = 0.1;
        solverObj = [];
        ic = [];
    end
    methods
        function obj = solver_system(varargin)
            warning("off", "SystemBlock:MATLABSystem:ParameterWithUnspecifiedDefaultValueRestrictedToBuiltinDataType"); 
            obj@matlab.System(varargin{:});
        end
    end
    methods(Access = protected)
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants
            if ~isempty(obj.ic)
                obj.solverObj.set_initial_states(obj.ic);
            end
        end

        function num = getNumOutputsImpl(obj)
            num = 1; % concentration
            % csd
            if ~isempty(obj.solverObj)
                num = num + numel(obj.solverObj.props);
            end
        end

        function [c, varargout] = getOutputNamesImpl(obj)
            % Return output port names for System block
            c = 'c';
            if ~isempty(obj.solverObj)
                n = numel(obj.solverObj.props);
                varargout = cell(n, 1);
                for i = 1 : n
                    varargout{i} = ['csd' num2str(i)];
                end
            end
        end

        function [c, varargout] = stepImpl(obj,u)
            % Implement algorithm. Calculate y as a function of input u and
            % discrete states.
            in = make_inputs(u);
            nextState = obj.solverObj.step([], in, obj.cfl, obj.sampleTimeSecond);
            c = nextState(1).conc;
            varargout = {nextState.csd};
        end

        function [c, varargout] = getOutputSizeImpl(obj)
            c = [1 1];
            if ~isempty(obj.solverObj)
                n = numel(obj.solverObj.props);
                varargout = cell(n, 1);
                for i = 1 : n
                    varargout{i} = [obj.solverObj.props(1).sizeGrids.nPoints 1];
                end
            end
        end

        function [c, varargout] = getOutputDataTypeImpl(obj)
            c = "double";
            if ~isempty(obj.solverObj)
                n = numel(obj.solverObj.props);
                sa = repmat("double", n, 1);
                varargout = sa.cellstr;
            end
        end

        function [c, varargout] = isOutputComplexImpl(obj)
            c = false;
            if ~isempty(obj.solverObj)
                n = numel(obj.solverObj.props);
                varargout = repmat({false}, n, 1);
            end
        end

        function [c, varargout] = isOutputFixedSizeImpl(obj)
            c = true;
            if ~isempty(obj.solverObj)
                n = numel(obj.solverObj.props);
                varargout = repmat({true}, n, 1);
            end
        end
    end
end
