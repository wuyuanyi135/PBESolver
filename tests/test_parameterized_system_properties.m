classdef test_parameterized_system_properties < matlab.unittest.TestCase

    
    methods(Test)
        function test_solubility(testCase)
            props = parameterized_system_properties();
            testCase.assertEqual(props.solubility(40), 0.0193, 'AbsTol', 1e-3);
        end
        
        function test_default_kinetics(testCase)
            props = parameterized_system_properties();
            
            conc = 0.028914;
            svar = struct();
            svar.tC = 40;
            svar.vf = 0.1;
            sol = props.solubility(40);
            svar.sigma = conc / sol - 1;
            
            [GD, Bp, Bs] = props.kinetics(svar);
            testCase.assertEqual(Bp, 25000000, 'RelTol', 1e-3);
            testCase.assertEqual(Bs, 5.3861e+08, 'RelTol', 1e-3);
            testCase.assertEqual(GD, 0.05, 'RelTol', 1e-3);
            
            % dissolution
            conc = 0.009638;
            svar.sigma = conc / sol - 1;
            [GD, ~, ~] = props.kinetics(svar);
            testCase.assertEqual(GD, -1.1, 'RelTol', 1e-3);
        end
    end
end

