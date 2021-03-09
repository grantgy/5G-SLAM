classdef Testclass < matlab.unittest.TestCase
  methods (Tests)
        function test_kf_prediction_symmetrix_P(testCase)
            %test_kf_prediction_symmetrix_P tests that kf_prediction works as intended
            
            %This is to check the covariance matrix after Kalman filter prediction step 
            % should be symmetric 
            x = roundn(5*rand(4,1),-3);
            P = roundn(5*abs(rand(4,4)),-3);
            P = (P+P')/2;
            F = roundn(5*rand(4,4),-3);
            Q = roundn(5*abs(rand(4,4)),-3);
            Q = (Q+Q')/2;            
            [~,actSolution_P] = kf_prediction(x,P,F,Q);
            testCase.verifyEqual(actSolution_P,actSolution_P');
        end
        
        function test_kf_prediction_dimension(testCase)
            % test_kf_prediction_dimension tests that kf_prediction works as intended
            
            % This is to check the size of state and covariance matrix outputed by Kalman 
            % filter prediction step meet the size requirement.
            x = roundn(5*rand(4,1),-3);
            P = roundn(5*abs(rand(4,4)),-3);
            P = (P+P')/2;
            F = roundn(5*rand(4,4),-3);
            Q = roundn(5*abs(rand(4,4)),-3);
            Q = (Q+Q')/2;
            
            expSolution_x_dim = [4,1];
            expSolution_P_dim = [4,4];
            
            [actSolution_x,actSolution_P] = kf_prediction(x,P,F,Q);

            testCase.verifySize (actSolution_P,expSolution_P_dim);
            testCase.verifySize (actSolution_x,expSolution_x_dim);
        end
        
        function test_kf_prediction_case(testCase)
            % test_kf_prediction_case is a case test to test that kf_prediction works as intended
            
            x = [1;2];
            P = [1,0.5;0.5,1];
            F = [1,2;0,1];
            Q = [1,0;0,1];
            
            expSolution_x = [5;2];
            expSolution_P = [8,2.5;2.5,2];
            
            [actSolution_x,actSolution_P] = kf_prediction(x,P,F,Q);
            testCase.verifyEqual(actSolution_x,expSolution_x);
            testCase.verifyEqual(actSolution_P,expSolution_P);
        end
  
        function test_Cali_AOA_increase(testCase)
            % test_Cali_AOA_increase tests that Cali.m works as intended if
            % AOA needs to be increased
            
            actSolution = Cali([0;-pi/2;0;0;0],[0;3*pi/2;0;0;0]);
            expSolution = [0;3*pi/2;0;0;0];
            testCase.verifyEqual(actSolution,expSolution);
        end
        
        function test_Cali_AOA_decrease(testCase)
            % test_Cali_AOA_decrease tests that Cali.m works as intended if
            % AOA need to be decreased
            actSolution = Cali([0;3*pi/4;0;0;0],[0;-pi/2;0;0;0]);
            expSolution = [0;-5*pi/4;0;0;0];
            testCase.verifyEqual(actSolution,expSolution);
        end
        
        function test_Cali_AOA_unchange(testCase)
            % test_Cali_AOA_unchange tests that Cali.m works as intended if
            % AOA do not need to be changed
            actSolution = Cali([0;pi/2;0;0;0],[0;-pi/2;0;0;0]);
            expSolution = [0;pi/2;0;0;0];
            testCase.verifyEqual(actSolution,expSolution);
        end
        
        function test_Cali_AOD_increase(testCase)
            % test_Cali_AOD_increase tests that Cali.m works as intended if
            % AOD need to be increased
            actSolution = Cali([0;0;0;-pi/2;0],[0;0;0;3*pi/2;0]);
            expSolution = [0;0;0;3*pi/2;0];
            testCase.verifyEqual(actSolution,expSolution);
        end
  end
end
