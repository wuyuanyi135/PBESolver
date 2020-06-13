function options = make_options()
options = struct();
options.earlyStopThreshold = 1e-5;

% Enable sub-CSD optimization
options.useSubCSD = true;
end

