function [pos,Cov] = newBern(measurement,R,x,type)

    x = [x(1:2,1);0;x(3:end,1)];
    pos_b = [0,0,40]';
    S = chol(R,'lower');
    w = sqrt(5)*[1 0 0 0 0 -1  0  0  0 0;
                 0 1 0 0 0  0 -1  0  0 0;
                 0 0 1 0 0  0  0 -1  0 0;
                 0 0 0 1 0  0  0  0 -1 0;
                 0 0 0 0 1  0  0  0  0 -1];
    pos = zeros(3,1);
    Cov= zeros(3);
    point_cub = S*w+measurement;
    
    if type ==2
        for j =1:10            
            TOA = point_cub(1,j); AOAaz = point_cub(4,j); AOAel = point_cub(5,j);
            Radius = TOA-x(5);
            R_xy = Radius*cos(AOAel);
            pos_0 = [x(1) + R_xy*cos(AOAaz + x(4)); x(2) + R_xy*sin(AOAaz + x(4)); x(3) + Radius*sin(AOAel)];
            x_v = cubature_point_estimation_VA(point_cub(:,j),x,pos_b,R,pos_0);
            pos = pos + x_v;
            Cov = Cov + x_v*x_v';
        end
        pos = pos./10;
        Cov = diag(diag(Cov./10 - pos*pos'));
    elseif type ==3
        for j =1:10         
            TOA = point_cub(1,j); AOAaz = point_cub(4,j); AOAel = point_cub(5,j);
            Radius = TOA-x(5);
            R_xy = Radius*cos(AOAel);          
            VA_geo = [x(1) + R_xy*cos(AOAaz + x(4)); x(2) + R_xy*sin(AOAaz + x(4)); x(3) + Radius*sin(AOAel)]; 
            u = (pos_b-VA_geo)/norm(pos_b-VA_geo);
            f = (pos_b+VA_geo)/2;
            pos_0 = VA_geo + (f-VA_geo)'*u/( (x(1:3,1) - VA_geo)'*u )* (x(1:3,1)-VA_geo);            
            x_s = cubature_point_estimation_SP(point_cub(:,j),x,pos_b,R,pos_0);
            pos = pos + x_s;
            Cov = Cov + x_s*x_s';
        end
        pos = pos./10;
        Cov = diag(diag(Cov./10 - pos*pos'));
    else
        error('Wrong input source');
    end
end

function [pos_v] = cubature_point_estimation_VA(point_cub,x,pos_b,R,pos_0)

    pos_v =pos_0;
    C = eye(size(R))/R;
    cost = inf;
    cost = [cost,(VA_function(pos_v,x,pos_b)-point_cub)'*C*(VA_function(pos_v,x,pos_b)-point_cub)];
    cost(1) =cost(2)+1;    
    iter =1;
    beta = 0.2;
    
    while cost(iter)-cost(iter+1) > 0 && iter < 100
        h = VA_function(pos_v,x,pos_b);
        if beta+.2 <= 1
            beta = beta + .2;
        end        
        J = (Cali(VA_function(pos_v+[0.001;0;0],x,pos_b),h) -  h);
        J =[J, (Cali(VA_function(pos_v+[0;0.001;0],x,pos_b),h) -  h)];
        J =[J, (Cali(VA_function(pos_v+[0;0;0.001],x,pos_b),h) -  h)];
        J = J./0.001;
        h=Cali(h,point_cub);
        error = h-point_cub;        
        B = -J*pos_v + error;
        BM_hat = (J'*C*J + J'*C'*J)\(-J'*C*B-J'*C'*B);
        pos_v = (1-beta)*pos_v + beta*BM_hat;
        cost = [cost,error'*C*error];
        iter = iter + 1;
    end
end

function [pos_s] = cubature_point_estimation_SP(point_cub,x,pos_b,R,pos_0)

    pos_s = pos_0;
    C = eye(size(R))/R;
    cost = inf;
    cost = [cost,(SP_function(pos_s,x,pos_b)-point_cub)'/R*(SP_function(pos_s,x,pos_b)-point_cub)];
    cost(1) =cost(2)+1;
    iter =1;
    beta =0.2;
    
    while cost(iter)-cost(iter+1) > 0 && iter < 100
        h = SP_function(pos_s,x,pos_b);
        
        if beta+.2 <= 1
            beta = beta + .2;
        end
        
        J = (Cali(SP_function(pos_s+[0.001;0;0],x,pos_b),h) -  h);
        J =[J, (Cali(SP_function(pos_s+[0;0.001;0],x,pos_b),h) -  h)];
        J =[J, (Cali(SP_function(pos_s+[0;0;0.001],x,pos_b),h) -  h)];
        J = J./0.001;       
        h=Cali(h,point_cub);
        error = h-point_cub;
        B = -J*pos_s + error;
        BM_hat = (J'*C*J + J'*C'*J)\(-J'*C*B-J'*C'*B);
        pos_s = (1-beta)*pos_s + beta*BM_hat; 
        cost = [cost,error'/R*error];
        iter = iter + 1;
    end    
end

function [z] = VA_function(pos_v,x,pos_b)

    u = (pos_b -pos_v) ./norm(pos_b - pos_v);
    f = (pos_b +pos_v) ./2;
    pos_s= pos_v +(f-pos_v)'*u/((x(1:3,1)-pos_v)'*u)*(x(1:3,1)-pos_v);
    z = zeros(5,1);    
    z(1,1) = norm(pos_v-x(1:3,1))+x(5,1);
    z(2,1) = atan2(pos_s(2,1),pos_s(1,1));
    z(3,1) = asin((pos_s(3,1)-pos_b(3,1))/norm(pos_s(1:3,1)-pos_b(1:3,1)));
    z(4,1) = atan2(pos_v(2,1)-x(2,1),pos_v(1,1)-x(1,1)) - x(4,1);
    z(5,1) = asin((pos_v(3,1)-x(3,1))/norm(x(1:3,1)-pos_v(1:3,1)));
    
    while z(2,1) < 0
        z(2,1) = z(2,1)+2*pi;
    end
    
    while z(4,1) < 0
        z(4,1) = z(4,1)+2*pi;
    end
end 

function [z] = SP_function(pos_s,x,pos_b)

    z = zeros(5,1);
    z(1,1) = norm(pos_s-x(1:3,1))+norm(pos_s-pos_b(1:3,1))+x(5,1);
    z(2,1) = atan2(pos_s(2,1),pos_s(1,1));
    z(3,1) = asin((pos_s(3,1)-pos_b(3,1))/norm(pos_s(1:3,1)-pos_b(1:3,1)));
    z(4,1) = atan2(pos_s(2,1)-x(2,1),pos_s(1,1)-x(1,1)) - x(4,1);
    z(5,1) = asin((pos_s(3,1)-x(3,1))/norm(x(1:3,1)-pos_s(1:3,1)));
    
    while z(2,1) < 0
        z(2,1) = z(2,1)+2*pi;
    end
    
    while z(4,1) < 0
        z(4,1) = z(4,1)+2*pi;
    end
end 
