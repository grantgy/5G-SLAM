function [z] = measurement_function(pos,x,pos_b,m,measurement)
    period = 2*pi;
    z = zeros(5,1);
    x = [x(1:2,1);0;x(3:end,1)];

    if m ==1  
        z(1,1) = norm(pos_b(1:3,1)-x(1:3,1))+x(5,1);
        z(2,1) = atan2(x(2,1),x(1,1));
        z(3,1) = asin((x(3,1)-pos_b(3,1))/norm(x(1:3,1)-pos_b(1:3,1)));
        z(4,1) = pi + atan2(x(2,1),x(1,1)) - x(4,1);
        z(5,1) = asin((pos_b(3,1)-x(3,1))/norm(x(1:3,1)-pos_b(1:3,1)));  
    elseif m ==2
        pos_v = pos;
        u = (pos_b -pos_v) ./norm(pos_b - pos_v);
        f = (pos_b +pos_v) ./2;
        pos_s= pos_v +(f-pos_v)'*u/((x(1:3,1)-pos_v)'*u)*(x(1:3,1)-pos_v);  
        z(1,1) = norm(pos_v-x(1:3,1))+x(5,1);
        z(2,1) = atan2(pos_s(2,1),pos_s(1,1));
        z(3,1) = asin((pos_s(3,1)-pos_b(3,1))/norm(pos_s(1:3,1)-pos_b(1:3,1)));
        z(4,1) = atan2(pos_v(2,1)-x(2,1),pos_v(1,1)-x(1,1)) - x(4,1);
        z(5,1) = asin((pos_v(3,1)-x(3,1))/norm(x(1:3,1)-pos_v(1:3,1)));
    elseif m ==3
        pos_s = pos;
        z(1,1) = norm(pos_s-x(1:3,1))+norm(pos_s-pos_b(1:3,1))+x(5,1);
        z(2,1) = atan2(pos_s(2,1),pos_s(1,1));
        z(3,1) = asin((pos_s(3,1)-pos_b(3,1))/norm(pos_s(1:3,1)-pos_b(1:3,1)));
        z(4,1) = atan2(pos_s(2,1)-x(2,1),pos_s(1,1)-x(1,1)) - x(4,1);
        z(5,1) = asin((pos_s(3,1)-x(3,1))/norm(x(1:3,1)-pos_s(1:3,1)));
    else
        error('Wrong input source');
    end

    while z(2,1) < 0
        z(2,1) = z(2,1)+period;
    end

    while z(4,1) < 0
        z(4,1) = z(4,1)+period;
    end

    while z(4,1) >= 2*pi
        z(4,1) = z(4,1)-period;
    end

    z = Cali(z,measurement);
end 
