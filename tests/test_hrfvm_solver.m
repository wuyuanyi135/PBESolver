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
    end
end

