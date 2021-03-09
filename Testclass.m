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
  
  end
end
