classdef test_size_grid < matlab.unittest.TestCase
   
    %% Test Method Block
    methods (Test)
        function test_to_array(testCase)
           s = size_grid(0, 999, 1000);
           testCase.assertEqual(s.to_array(), (0:999)'); 
           
        end
        
        function test_interval(testCase)
           s = size_grid(0, 999, 1000);
           testCase.assertEqual(s.interval(), 1); 
           
        end
    end
end