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
    end
end