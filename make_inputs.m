function inputs = make_inputs(tC, resTime, inConc, inCSDs)
arguments
    tC = 25;
    resTime = inf;
    inConc = 0;
    inCSDs = {};
end
inputs = struct();
inputs.tC = tC;
inputs.resTime = resTime;
inputs.inConc = inConc;
inputs.inCSDs = inCSDs;
end

