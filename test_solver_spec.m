classdef test_solver_spec < matlab.unittest.TestCase
    % Unit test to ensure the solver satisfies the solver fundamental
    % specification or assumption.
    
    properties(TestParameter)
        useSubCSD = struct('UseSubCSD', true, 'NoSubCSD', false);
        solver = struct('Upwind', @upwind_solver, 'HRFVM', @hrfvm_solver);
    end
    
    methods(Test)
        function test_growth_touch_boundary(testCase, useSubCSD, solver)
            % When the solving steps touch the size grid boundary (right),
            % the boundary should not affect the result of the other
            % points.
            % If this test passes, the boundary of the CSD, or the boundary
            % of the subCSD does not affect those channel adjecent to the
            % boundary. Therefore, cutting at the effective size should
            % only affect the cut channels but does not affect those
            % channel smaller than the cut.
            props = parameterized_system_properties();
            props.pnKp = 0;
            props.snKb = 0;
            g = size_grid(99, 100);
            props.sizeGrids = g;
            options = make_options();
            options.useSubCSD = useSubCSD;
            
            s = solver(options, props);
            
            nonBoundaryCenter = 50;
            nonBoundaryIc = make_states(props.solubility(60), normpdf(g.to_array(), nonBoundaryCenter, 5), props);
            boundaryCenter = 90;
            boundaryIc = make_states(props.solubility(60), normpdf(g.to_array(), boundaryCenter, 5), props);

            in = make_inputs(25);
            
            for tSpan = 25:25:100
                nonBoundaryNextState = s.step(nonBoundaryIc, in, tSpan);
                
                boundaryNextState = s.step(boundaryIc, in, tSpan);
                
                testCase.assertEqual( ...
                    nonBoundaryNextState.csd(nonBoundaryCenter-5:nonBoundaryCenter+5), ...
                    boundaryNextState.csd(boundaryCenter-5:boundaryCenter+5), ...
                    'RelTol', 1e-3 ...
                    );
            end
        end
        function test_dissolution_touch_boundary(testCase, useSubCSD, solver)
            % Refer to test_growth_touch_boundary.
            props = parameterized_system_properties();
            props.pnKp = 0;
            props.snKb = 0;
            g = size_grid(99, 100);
            props.sizeGrids = g;
            options = make_options();
            options.useSubCSD = useSubCSD;
            
            s = solver(options, props);
            
            nonBoundaryCenter = 50;
            nonBoundaryIc = make_states(props.solubility(25), normpdf(g.to_array(), nonBoundaryCenter, 5), props);
            boundaryCenter = 10;
            boundaryIc = make_states(props.solubility(25), normpdf(g.to_array(), boundaryCenter, 5), props);

            in = make_inputs(40);
            
            tSpan = 10;
            nonBoundaryNextState = s.step(nonBoundaryIc, in, tSpan);
            
            boundaryNextState = s.step(boundaryIc, in, tSpan);
            
            testCase.assertEqual( ...
                nonBoundaryNextState.csd(nonBoundaryCenter-5:nonBoundaryCenter+5), ...
                boundaryNextState.csd(boundaryCenter-5:boundaryCenter+5), ...
                'RelTol', 1e-3 ...
                );
            
        end
        
        
        function test_growth_only(testCase, useSubCSD, solver)
            props = parameterized_system_properties;
            props.pnKp = 0;
            props.snKb = 0;
            options = make_options();
            options.useSubCSD = useSubCSD;
            s = solver(options, props);
            sizeGrid = props.sizeGrids;
            g = sizeGrid.to_array();
            ic = make_states(props.solubility(60), normpdf(g, 20, 10), props);
            in = make_inputs(25);
            tSpan = 300;
            nextState = s.step(ic, in, tSpan);
            
            % Peak position
            [~, I] = max(nextState.csd);
            testCase.assertEqual(I, 95, 'AbsTol', 1);
        end
        
        function test_dissolution_only(testCase, useSubCSD, solver)
            props = parameterized_system_properties;
            props.pnKp = 0;
            props.snKb = 0;
            options = make_options();
            options.useSubCSD = useSubCSD;
            s = solver(options, props);
            sizeGrid = props.sizeGrids;
            g = sizeGrid.to_array();
            ic = make_states(props.solubility(25), normpdf(g, 50, 10), props);
            in = make_inputs(60);
            tSpan = 10;
            nextState = s.step(ic, in, tSpan);
            
            % Peak position
            [~, I] = max(nextState.csd);
            testCase.assertEqual(I, 35, 'AbsTol', 1);
        end
        
        function test_sub_csd_different_step(testCase, solver)
            % This test checks with different time step scales, whether the
            % concentration profile (by comparing the endpoints) are the
            % same. Practically, it verifies that under the circumstance of 
            % enabled sub-CSD optimization, whether the problem of leaking
            % channel is mitigated so that the result is reproducible at
            % different time step scale.
            scales = [0.3, 0.1, 0.05, 0.01];
            cs = zeros('like', scales);
            for i = 1: numel(scales)
                props = parameterized_system_properties;
                options = make_options();
                options.useSubCSD = 0; % should be true
                options.timeStepScale = scales(i);
                s = solver(options, props);
                ic = make_states(props.solubility(60), 0, props);
                in = make_inputs(25);
                tSpan = 500;
                nextState = s.step(ic, in, tSpan);
                cs(i) = nextState.conc;
            end
            
            stdDev = std(cs);
            mu = repmat(mean(cs), size(cs));
            
            testCase.assertEqual(mu, cs, 'AbsTol', 3*stdDev);
        end
    end
end

