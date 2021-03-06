function [pos,P,z_inv,R_inv] = CKF(x,measurement,pos,P,R,type)
    pos_b = [0,0,40]';
    x = [x(1:2,1);0;x(3:end,1)];
    S = chol(P,'lower');
    
    w = sqrt(size(pos_b,1))*[eye(size(pos_b,1)),-eye(size(pos_b,1))];
    
    pos_u = S*w+pos;
    z = zeros(5,6);
    z_inv = zeros(5,1);
    R_inv = zeros(5);
    P_xz  = zeros(3,5);
    if type == 1
        for i = 1:size(w,2)
            z(1,i) = norm(pos_u(1:3,i)-x(1:3,1))+x(5,1);
            z(2,i) = atan2(x(2,1),x(1,1));
            z(3,i) = asin((x(3,1)-pos_u(3,i))/norm(x(1:3,1)-pos_u(1:3,i)));
            z(4,i) = pi + atan2(x(2,1),x(1,1)) - x(4,1);
            z(5,i) = asin((pos_u(3,i)-x(3,1))/norm(x(1:3,1)-pos_u(1:3,i)));
            
            
            z(:,i) = calibration(z(:,i),z(:,1));
            
            z_inv = z_inv + z(:,i);
            R_inv = R_inv + z(:,i)*z(:,i)' ;
        end
        z_inv = z_inv./size(w,2);
        R_inv = R_inv./size(w,2) - z_inv*z_inv' + R;
        z_inv = Cali(z_inv, measurement);
        
    elseif type == 2 
        for i = 1:size(w,2)
            u = (pos_b -pos_u(:,i)) ./norm(pos_b - pos_u(:,i));
            f = (pos_b +pos_u(:,i)) ./2;
            pos_s= pos_u(:,i) +(f-pos_u(:,i))'*u/((x(1:3,1)-pos_u(:,i))'*u)*(x(1:3,1)-pos_u(:,i));
            z(1,i) = norm(pos_u(1:3,i)-x(1:3,1))+x(5,1);
            z(2,i) = atan2(pos_s(2,1),pos_s(1,1));
            z(3,i) = asin((pos_s(3,1)-pos_b(3,1))/norm(pos_s(1:3,1)-pos_b(1:3,1)));
            z(4,i) = atan2(pos_u(2,i)-x(2,1),pos_u(1,i)-x(1,1)) - x(4,1);
            z(5,i) = asin((pos_u(3,i)-x(3,1))/norm(x(1:3,1)-pos_u(1:3,i)));     
            
            z(:,i) = calibration(z(:,i),z(:,1));
            
            z_inv = z_inv + z(:,i);
            R_inv = R_inv + z(:,i)*z(:,i)';
            P_xz  = P_xz  + pos_u(:,i)*z(:,i)';
        end
        [pos,P,z_inv,R_inv]= update_parameters (pos,P,z_inv,R_inv,P_xz,R,measurement);

        
    elseif type == 3 
        for i = 1:size(w,2)
            z(1,i) = norm(pos_u(1:3,i)-x(1:3,1))+norm(pos_u(1:3,i)-pos_b(1:3,1))+x(5,1);
            z(2,i) = atan2(pos_u(2,i),pos_u(1,i));
            z(3,i) = asin((pos_u(3,i)-pos_b(3,1))/norm(pos_u(1:3,i)-pos_b(1:3,1)));
            z(4,i) = atan2(pos_u(2,i)-x(2,1),pos_u(1,i)-x(1,1)) - x(4,1);
            z(5,i) = asin((pos_u(3,i)-x(3,1))/norm(x(1:3,1)-pos_u(1:3,i)));       
            
            z(:,i) = calibration(z(:,i),z(:,1));
            
            z_inv = z_inv + z(:,i);
            R_inv = R_inv + z(:,i)*z(:,i)';
            P_xz  = P_xz  + pos_u(:,i)*z(:,i)';
        end
        [pos,P,z_inv,R_inv]= update_parameters (pos,P,z_inv,R_inv,P_xz,R,measurement);
        
    else
        error('Wrong input source');
    end
    
end

function output = calibration (input,reference)
    while input(2,1) < 0
        input(2,1) = input(2,1)+2*pi;
    end
            
    while input(4,1) < 0
        input(4,1) = input(4,1)+2*pi;
    end
	output = Cali(input,reference);
end


function [pos,P,z_inv,R_inv]= update_parameters (pos,P,z_inv,R_inv,P_xz,R,measurement)
    z_inv = z_inv./6;
    R_inv = R_inv./6 - z_inv*z_inv' + R;
    P_xz  = P_xz./6  - pos*z_inv';
    K     = P_xz  / R_inv;
    z_inv = Cali(z_inv, measurement);
    pos   = pos   + K*(measurement - z_inv);
    P     = P     - K*R_inv*K';
end
