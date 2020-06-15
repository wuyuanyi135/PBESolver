classdef test_size_grid < matlab.unittest.TestCase
   
    %% Test Method Block
    methods (Test)
        function test_to_array(testCase)
           s = size_grid(999, 1000);
           testCase.assertEqual(s.to_array(), (0:999)'); 
           
        end
        
        function test_interval(testCase)
           s = size_grid(999, 1000);
           testCase.assertEqual(s.interval(), 1); 
        end
        
        function test_equality(testCase)
            s1 = size_grid(999, 1000);
            s2 = size_grid(999, 1000);
            testCase.assertEqual(s1, s2);
            
            s3 = size_grid(992, 1000);
            testCase.assertNotEqual(s1, s3);
            
            s4 = size_grid(999, 999);
            testCase.assertNotEqual(s1, s4);
        end
    end
end