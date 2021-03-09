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
        
        function test_Cali_AOD_decrease(testCase)
            % test_Cali_AOD_decrease tests that Cali.m works as intended if
            % AOD need to be decreased
            actSolution = Cali([0;0;0;3*pi/4;0],[0;0;0;-pi/2;0]);
            expSolution = [0;0;0;-5*pi/4;0];
            testCase.verifyEqual(actSolution,expSolution);
        end
        
        function test_Cali_AOD_unchange(testCase)
            % test_Cali_AOD_unchange tests that Cali.m works as intended if
            % AOD do not need to be changed
            actSolution = Cali([0;pi/3;0;0;0],[0;-pi/2;0;0;0]);
            expSolution = [0;pi/3;0;0;0];
            testCase.verifyEqual(actSolution,expSolution);
        end
        
        function test_UM(testCase)
            % test_UM tests UM.m can generate correct class of uniform
            % distribution mixture
            expSolution_weight = [0.5,0.5];
            expSolution_pro = [4,3];            
            actSolution_UM = UM([0.5,0.5],[4,3]);
            testCase.verifyEqual(actSolution_UM.weight,expSolution_weight);
            testCase.verifyEqual(actSolution_UM.pro,expSolution_pro);
        end
        
        function test_GM(testCase)
            % test_GM tests GM.m can generate correct class of Gaussian
            % mixture
            expSolution_weight = 1;
            expSolution_mean = [0;0;40]; 
            expSolution_covariance = eye(3);
            actSolution_GM = GM(1,[0;0;40],eye(3));
            testCase.verifyEqual(actSolution_GM.weight,expSolution_weight);
            testCase.verifyEqual(actSolution_GM.mean,expSolution_mean);
            testCase.verifyEqual(actSolution_GM.covariance,expSolution_covariance);
        end
        
        function test_MAP(testCase)
            % test_MAP tests MAP.m can generate correct class of map of
            % landmarks
            expSolution_BS = [0;0;0];
            expSolution_VA = [200,100;200,100;40,20]; 
            expSolution_SP = [50,-50,30;50,-50,30;0,10,20];
            actSolution_MAP = MAP([0;0;0],[200,100;200,100;40,20],[50,-50,30;50,-50,30;0,10,20]);
            testCase.verifyEqual(actSolution_MAP.BS,expSolution_BS);
            testCase.verifyEqual(actSolution_MAP.VA,expSolution_VA);
            testCase.verifyEqual(actSolution_MAP.SP,expSolution_SP);
        end
        
        function test_MB(testCase)
            % test_MB tests MB.m can generate correct class of
            % Multi-Bernoulli process
            expSolution_p_exist = [1,0.5];
            expSolution_Gaussian_mixture = [GM(1,[0;0;40],eye(3)),GM(0.1,[0;0;20],eye(3))]; 
            expSolution_log_weight = [100,-2.1];
            expSolution_birth_time = 0;
            
            actSolution_MB = MB([1,0.5],[GM(1,[0;0;40],eye(3)),GM(0.1,[0;0;20],eye(3))],[100,-2.1],0);
            testCase.verifyEqual(actSolution_MB.p_exist,expSolution_p_exist);
            testCase.verifyEqual(actSolution_MB.Gaussian_mixture,expSolution_Gaussian_mixture);
            testCase.verifyEqual(actSolution_MB.log_weight,expSolution_log_weight);
            testCase.verifyEqual(actSolution_MB.birth_time,expSolution_birth_time);
        end
  end
end
