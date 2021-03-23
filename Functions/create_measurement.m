function [measurement,v_state] = create_measurement(x,steps,v,ang_v,T,R,L_BS,L_VA,L_SP)
    x = [x(1:2,1);0;x(3:end,1)];

    measurement = cell(1,steps);
    P_detection = 0.9; % detection probability
    v_state = x;
    
    for i =1:steps
       v_a = v/ang_v;
       ang_delta = ang_v * T;
        
       x(1,1) = x(1,1) + v_a*(sin(x(4,1)+ang_delta) -sin(x(4,1)));
       x(2,1) = x(2,1) + v_a*(-cos(x(4,1)+ang_delta) +cos(x(4,1)));
       x(4,1) = x(4,1) + ang_delta;
        
       v_state = [v_state,x];
           
       z1 = [];
       for j=1:size(L_BS,2)
           
           if rand(1)<P_detection
               z = zeros(5,1);
               z(1,1) = norm(L_BS(1:3,j)-x(1:3,1))+x(5,1);
               z(2,1) = atan2(x(2,1),x(1,1));
               z(3,1) = asin((x(3,1)-L_BS(3,j))/norm(x(1:3,1)-L_BS(1:3,j)));
               z(4,1) = pi + atan2(x(2,1),x(1,1)) - x(4,1);
               z(5,1) = asin((L_BS(3,j)-x(3,1))/norm(x(1:3,1)-L_BS(1:3,j)));    
               
               z1 = [z1,z+mvnrnd([0;0;0;0;0],R,1)'];
           end 
       end
       
       for j=1:size(L_VA,2)
           if rand(1)<P_detection
               z = zeros(5,1);
               u = (L_BS -L_VA(:,j)) ./norm(L_BS - L_VA(:,j));
               f = (L_BS +L_VA(:,j)) ./2;
               pos_s= L_VA(:,j) +(f-L_VA(:,j))'*u/((x(1:3,1)-L_VA(:,j))'*u)*(x(1:3,1)-L_VA(:,j));
               z(1,1) = norm(L_VA(1:3,j)-x(1:3,1))+x(5,1);
               z(2,1) = atan2(pos_s(2,1),pos_s(1,1));
               z(3,1) = asin((pos_s(3,1)-L_BS(3,1))/norm(pos_s(1:3,1)-L_BS(1:3,1)));
               z(4,1) = atan2(L_VA(2,j)-x(2,1),L_VA(1,j)-x(1,1)) - x(4,1);
               z(5,1) = asin((L_VA(3,j)-x(3,1))/norm(x(1:3,1)-L_VA(1:3,j)));
               
               z1 = [z1,z+mvnrnd([0;0;0;0;0],R,1)'];
           end
       end
       
       for j=1:size(L_SP,2)
           if rand(1)<P_detection
               if norm(L_SP(1:3,j)-x(1:3,1)) <= 50
                   z = zeros(5,1);
                   z(1,1) = norm(L_SP(1:3,j)-x(1:3,1))+norm(L_SP(1:3,j)-L_BS(1:3,1))+x(5,1);
                   z(2,1) = atan2(L_SP(2,j),L_SP(1,j));
                   z(3,1) = asin((L_SP(3,j)-L_BS(3,1))/norm(L_SP(1:3,j)-L_BS(1:3,1)));
                   z(4,1) = atan2(L_SP(2,j)-x(2,1),L_SP(1,j)-x(1,1)) - x(4,1);
                   z(5,1) = asin((L_SP(3,j)-x(3,1))/norm(x(1:3,1)-L_SP(1:3,j)));
                   z1 = [z1,z+mvnrnd([0;0;0;0;0],R,1)'];
               end
           end 
       end
       for k = 1 :size(z1,2)
           while z1(2,k) < 0
               z1(2,k) = z1(2,k)+2*pi;
           end
           
           while z1(4,k) < 0
               z1(4,k) = z1(4,k)+2*pi;
           end
       end
       
       measurement{1,i} = z1;
    end
    v_state(3,:)=[]; 
end
