classdef test_hrfvm_solver < matlab.unittest.TestCase
    methods(Test)
        function test_default_growth(testCase)
            props = parameterized_system_properties();
            options = make_options();
            options.earlyStopThreshold = 0;
            s = hrfvm_solver(options, props);
            ic = make_states(props.solubility(60), 0, props);
            in = make_inputs(25);
            cfl = 0.1;
            tSpan = 300;
            nextState = s.step(ic, in, cfl, tSpan);
            testCase.assertEqual(nextState.conc, 0.030478231269330, 'RelTol', 1e-3);
        end
        
        function test_multistep(testCase)
            props = parameterized_system_properties();
            options = make_options();
            options.earlyStopThreshold = 0;
            s = hrfvm_solver(options, props);
            ic = make_states(props.solubility(60), 0, props);
            in = make_inputs(25);
            cfl = 0.1;
            tSpan = 100;
            nextState = s.step(ic, in, cfl, tSpan);
            testCase.assertEqual(nextState.conc, 0.036689862918977, 'AbsTol', 1e-4);

            nextState = s.step(nextState, in, cfl, tSpan);
            testCase.assertEqual(nextState.conc, 0.035507179592030, 'AbsTol', 1e-4);

            nextState = s.step(nextState, in, cfl, tSpan);
            testCase.assertEqual(nextState.conc, 0.030478231269330, 'AbsTol', 1e-4);
        end
        
        function test_early_stop(testCase)
            props = parameterized_system_properties();
            options = make_options();
            options.earlyStopThreshold = 1e-5;
            s = hrfvm_solver(options, props);
            ic = make_states(props.solubility(60), 0, props);
            in = make_inputs(25);
            cfl = 0.1;
            tSpan = 2000;
            nextState = s.step(ic, in, cfl, tSpan);
            testCase.assertEqual(nextState.conc, 0.010595206881399, 'AbsTol', 1e-4);
        end
        
        function test_polymorphism(testCase)
            props1 = parameterized_system_properties();
            props2 = parameterized_system_properties();
            props2.solubilityPoly = [6.222e-3; -1.165e-4; 7.644e-6];
            props = [props1, props2]; % TODO: props must be col vec!
            
            options = make_options();
            options.earlyStopThreshold = 0;
            
            ic = [
                make_states(props1.solubility(60), 0, props1),...
                make_states(nan, 0, props2)
            ];
        
            in = make_inputs(25);
            
            s = hrfvm_solver(options, props);
            
            cfl = 0.1;
            tSpan = 2000;
            nextState = s.step(ic, in, cfl, tSpan);            
        end
        
        function test_stateful_single_step(testCase)
            props = parameterized_system_properties();
            options = make_options();
            options.earlyStopThreshold = 0;
            s = hrfvm_solver(options, props);
            ic = make_states(props.solubility(60), 0, props);
            in = make_inputs(25);
            cfl = 0.1;
            tSpan = 300;
            s.set_initial_states(ic);
            nextState = s.step([], in, cfl, tSpan);
            testCase.assertEqual(nextState.conc, 0.030478231269330, 'RelTol', 1e-3);
            testCase.assertEqual(s.prevStates.conc, 0.030478231269330, 'RelTol', 1e-3);
        end
    end
end

