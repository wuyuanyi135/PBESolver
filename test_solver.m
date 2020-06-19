classdef test_solver < matlab.perftest.TestCase
    properties(TestParameter)
        useSubCSD = struct('UseSubCSD', true, 'NoSubCSD', false);
        solver = struct('Upwind', @upwind_solver, 'HRFVM', @hrfvm_solver);
    end
    methods(Test)
        function test_default_growth(testCase, useSubCSD, solver)
            props = parameterized_system_properties();
            options = make_options();
            options.earlyStopThreshold = 0;
            options.useSubCSD = useSubCSD;
            options.timeStepScale = 0.1;
            s = solver(options, props);
            ic = make_states(props.solubility(60), 0, props);
            in = make_inputs(25);
            tSpan = 300;
            nextState = s.step(ic, in, tSpan);
            % consistency test
            % Note, different solver, different options will result in
            % different final concentration. This is only for testing the
            % consistency but not ensuring the accuracy.
            % Note: Values hard-coded on 6/19/2020
            solverFcn = functions(solver);
            switch (solverFcn.function)
                case 'hrfvm_solver'
                    switch (useSubCSD)
                        case 1
                            testCase.assertEqual(nextState.conc, 0.0270545060780174, 'RelTol', 1e-3);
                        case 0
                            testCase.assertEqual(nextState.conc, 0.02649, 'RelTol', 1e-3);
                    end
                case 'upwind_solver'
                    switch (useSubCSD)
                        case 1
                            testCase.assertEqual(nextState.conc, 0.027965276284955, 'RelTol', 1e-3);
                        case 0
                            testCase.assertEqual(nextState.conc, 0.025437263363723, 'RelTol', 1e-3);
                    end
                otherwise
                    msg = sprintf('Solver %s does not have assertion case', solverFcn.function);
                    testCase.assertFail(msg);
            end
        end
        
        function test_multistep(testCase, useSubCSD, solver)
            % Note, this test uses higher allowed relative error because
            % multistep seems introduces accumulating error.
            props = parameterized_system_properties();
            options = make_options();
            options.earlyStopThreshold = 0;
            options.useSubCSD = useSubCSD;
            s = solver(options, props);
            ic = make_states(props.solubility(60), 0, props);
            in = make_inputs(25);
            tSpan = 100;
            nextState = s.step(ic, in, tSpan);
            nextState = s.step(nextState, in, tSpan);
            nextState = s.step(nextState, in, tSpan);
            
            solverFcn = functions(solver);
            switch (solverFcn.function)
                case 'hrfvm_solver'
                    switch (useSubCSD)
                        case 1
                            testCase.assertEqual(nextState.conc, 0.0270545060780174, 'RelTol', 1e-2);
                        case 0
                            testCase.assertEqual(nextState.conc, 0.02649, 'RelTol', 1e-2);
                    end
                case 'upwind_solver'
                    switch (useSubCSD)
                        case 1
                            testCase.assertEqual(nextState.conc, 0.027965276284955, 'RelTol', 1e-2);
                        case 0
                            testCase.assertEqual(nextState.conc, 0.025437263363723, 'RelTol', 1e-2);
                    end
                otherwise
                    msg = sprintf('Solver %s does not have assertion case', solverFcn.function);
                    testCase.assertFail(msg);
            end
        end
        
        function test_early_stop(testCase, useSubCSD, solver)
            props = parameterized_system_properties();
            options = make_options();
            options.earlyStopThreshold = 1e-5;
            options.useSubCSD = useSubCSD;
            s = solver(options, props);
            ic = make_states(props.solubility(60), 0, props);
            in = make_inputs(25);
            tSpan = 2000;
            nextState = s.step(ic, in, tSpan);
            testCase.assertEqual(nextState.conc, 0.01059, 'RelTol', 1e-3);
        end
        
        function test_polymorphism(testCase, useSubCSD, solver)
            props1 = parameterized_system_properties();
            props2 = parameterized_system_properties();
            props2.solubilityPoly = [6.222e-3; -1.165e-4; 7.644e-6];
            props = [props1, props2]; % TODO: props must be col vec!
            
            options = make_options();
            options.useSubCSD = useSubCSD;
            options.earlyStopThreshold = 0;
            
            ic = [
                make_states(props1.solubility(60), 0, props1),...
                make_states(nan, 0, props2)
            ];
        
            in = make_inputs(25);
            
            s = solver(options, props);
            tSpan = 2000;
            nextState = s.step(ic, in, tSpan); 
            testCase.assertEqual(nextState(1).conc, 0.008087, 'RelTol', 1e-3);
        end
        
        function test_stateful_single_step(testCase, useSubCSD, solver)
            props = parameterized_system_properties();
            options = make_options();
            options.earlyStopThreshold = 0;
            options.useSubCSD = useSubCSD;
            s = solver(options, props);
            ic = make_states(props.solubility(60), 0, props);
            in = make_inputs(25);
            tSpan = 300;
            s.set_initial_states(ic);
            nextState = s.step([], in, tSpan);
            conc = nextState.conc;
            prevConc = s.prevStates.conc;
            testCase.assertEqual(conc, prevConc);
            
            solverFcn = functions(solver);
            switch (solverFcn.function)
                case 'hrfvm_solver'
                    switch (useSubCSD)
                        case 1
                            testCase.assertEqual(nextState.conc, 0.0270545060780174, 'RelTol', 1e-3);
                        case 0
                            testCase.assertEqual(nextState.conc, 0.02649, 'RelTol', 1e-3);
                    end
                case 'upwind_solver'
                    switch (useSubCSD)
                        case 1
                            testCase.assertEqual(nextState.conc, 0.027965276284955, 'RelTol', 1e-3);
                        case 0
                            testCase.assertEqual(nextState.conc, 0.025437263363723, 'RelTol', 1e-3);
                    end
                otherwise
                    msg = sprintf('Solver %s does not have assertion case', solverFcn.function);
                    testCase.assertFail(msg);
            end
        end
        
        function test_size_dependent_growth(testCase, useSubCSD, solver)
            props = parameterized_system_properties();
            props.sizeGrids = size_grid(299, 300);
            props.gBeta = -2e-3;
            options = make_options();
            options.earlyStopThreshold = 0;
            options.useSubCSD = useSubCSD;
            s = solver(options, props);
            ic = make_states(props.solubility(60), 0, props);
            in = make_inputs(25);
            tSpan = 300;
            nextState = s.step(ic, in, tSpan);
            solverFcn = functions(solver);
            switch (solverFcn.function)
                case 'hrfvm_solver'
                    switch (useSubCSD)
                        case 1
                            testCase.assertEqual(nextState.conc, 0.032541771168378, 'RelTol', 1e-3);
                        case 0
                            testCase.assertEqual(nextState.conc, 0.032447564014222, 'RelTol', 1e-3);
                    end
                case 'upwind_solver'
                    switch (useSubCSD)
                        case 1
                            testCase.assertEqual(nextState.conc, 0.0289843825219405, 'RelTol', 1e-3);
                        case 0
                            testCase.assertEqual(nextState.conc, 0.026680361271864, 'RelTol', 1e-3);
                    end
                otherwise
                    msg = sprintf('Solver %s does not have assertion case', solverFcn.function);
                    testCase.assertFail(msg);
            end
        end
        function test_size_dependent_growth_increasing_faster(testCase, useSubCSD, solver)
            % Warning: this case needs attention. See TODO
            % size_dependent_mismatch.
            props = parameterized_system_properties();
            props.gBeta = 1e-3;
            options = make_options();
            options.earlyStopThreshold = 0;
            options.useSubCSD = useSubCSD;
            options.timeStepScale = 1e-1;
            s = solver(options, props);
            ic = make_states(props.solubility(60), 0, props);
            in = make_inputs(25);
            tSpan = 300;
            nextState = s.step(ic, in, tSpan);
            solverFcn = functions(solver);
            switch (solverFcn.function)
                case 'hrfvm_solver'
                    switch (useSubCSD)
                        case 1
                            testCase.assertEqual(nextState.conc, 0.031472683985705, 'RelTol', 1e-3);
                        case 0
                            testCase.assertEqual(nextState.conc, 0.0329554700906034, 'RelTol', 1e-3);
                    end
                case 'upwind_solver'
                    switch (useSubCSD)
                        case 1
                            testCase.assertEqual(nextState.conc, 0.027425509727680, 'RelTol', 1e-3);
                        case 0
                            testCase.assertEqual(nextState.conc, 0.024804926193285, 'RelTol', 1e-3);
                    end
                otherwise
                    msg = sprintf('Solver %s does not have assertion case', solverFcn.function);
                    testCase.assertFail(msg);
            end
        end
        function test_disallow_bad_growth_beta(testCase)
            props = parameterized_system_properties();
            props.sizeGrids = size_grid(999,1000);
            props.gBeta = -1e-3;
            function assign_beta()
                props.gBeta = -1e-2;
            end
            testCase.verifyError(@assign_beta, '');   
        end
    end
end

