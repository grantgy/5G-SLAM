function [particle,weight] = PF_prediction(particle,weight,v,ang_v,T,Q)
    v_a = v/ang_v;
    ang_delta = ang_v * T;
    

    for i = 1:size(particle,2)
        particle(1,i) = particle(1,i) + v_a*(sin(particle(3,i)+ang_delta) -sin(particle(3,i)));
        particle(2,i) = particle(2,i) + v_a*(-cos(particle(3,i)+ang_delta) +cos(particle(3,i)));
        particle(3,i) = particle(3,i) + ang_delta;
    end
    
    if size(particle,2) >1  % if there is only one particle, that is for test.
        particle = particle +  mvnrnd([0;0;0;0],Q,size(particle,2))';
    end
end