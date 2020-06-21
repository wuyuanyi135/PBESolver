function options = make_options()
options = struct();
options.earlyStopThreshold = 1e-5;

% Enable sub-CSD optimization
options.useSubCSD = true;

% computed maxmium time step's scale to the actual time step. (was cfl)
options.timeStepScale = 0.1;

% Is a MSMPR?
options.isMSMPR = false;

% Only works when isMSMPR is true.
options.residenceTimeStepScale = 0.1;
end

