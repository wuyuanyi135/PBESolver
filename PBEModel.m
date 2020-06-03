classdef PBEModel
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        solubility_fn
        kinetics_fn
        flux_limiter_fcn
        grad_ratio_fcn
        solver_fcn
    end
    
    methods
        function obj = PBEModel(solubility_fn, kinetics_fn, flux_limiter_fcn, grad_ratio_fcn, solver_fcn)
            obj.flux_limiter_fcn = flux_limiter_fcn;
            obj.solubility_fn = solubility_fn;
            obj.kinetics_fn = kinetics_fn;
            obj.grad_ratio_fcn = grad_ratio_fcn;
            obj.solver_fcn = solver_fcn;
        end
        
    end
end

