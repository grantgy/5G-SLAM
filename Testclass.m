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
  
  end
end
