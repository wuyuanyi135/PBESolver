classdef size_grid
    % size grid of a solution
    
    properties
        lowBound
        highBound
        nPoints
    end
    
    methods
        function obj = size_grid(highBound, nPoints)
            obj.lowBound = 0;
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
        
        function equality = eq(obj1, obj2)
            equality = obj1.highBound == obj2.highBound && obj1.interval == obj2.interval;
        end
    end
end

