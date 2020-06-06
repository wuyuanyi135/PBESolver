classdef size_grid
    % size grid of a solution
    
    properties
        lowBound
        highBound
        nPoints
    end
    
    methods
        function obj = size_grid(lowBound, highBound, nPoints)
            obj.lowBound = lowBound;
            obj.highBound = highBound;
            obj.nPoints = nPoints;
        end
        
        function a = to_array(obj)
            % Get the column vector format of the grid.
            % numel(obj) must equal 1.
            a = linspace(obj.lowBound, obj.highBound, obj.nPoints)';
        end
        
        function ret = interval(obj)
            % Get the grid interval in um
            ret = ([obj.highBound] - [obj.lowBound]) ./ ([obj.nPoints] - 1);
        end
    end
end

