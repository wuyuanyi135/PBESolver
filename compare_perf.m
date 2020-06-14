function [h] = compare_perf(baseline, new, varargin)
    filteredBaseline = [];
    filteredNew = [];
    
    targetMap = containers.Map({new.Name}, 1:numel(new));
    
    for i = 1 : numel(baseline)
        baseName = baseline(i).Name;
        if (targetMap.isKey(baseName))
            filteredNew = [filteredNew; new(targetMap(baseName))];
            filteredBaseline = [filteredBaseline; baseline(i)];
        end
    end
    
    h = comparisonPlot(filteredBaseline, filteredNew, varargin{:});
end

