classdef test_make_states < matlab.unittest.TestCase
   
    %% Test Method Block
    methods (Test)
        function test_empty(testCase)
            props = parameterized_system_properties;
            s = make_states(0.1, 0, props);
            testCase.assertEqual(s.conc, 0.1);
            testCase.assertEqual(s.csd, zeros(size(props.sizeGrids.to_array())));
            testCase.assertEqual(s.moment3, 0);
        end
        
        function test_uniform_seed(testCase)
            props = parameterized_system_properties;
            g = props.sizeGrids.to_array();
            seed = normpdf(g, 10, 10);
            s = make_states(0.1, seed, props);
            testCase.assertEqual(s.moment3, 4.091291360282422e+03, 'RelTol', 1e-3); 
        end
    end
end