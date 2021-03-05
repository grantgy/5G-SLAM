function [output] = Cali(output,input)
    
    while (output(2) - input(2) > pi)
        output(2) = output(2) - 2*pi;
    end
    
    while (output(2) - input(2) < -pi)
        output(2) = output(2) + 2*pi;
    end
    
    while (output(4) - input(4) > pi)
        output(4) = output(4) - 2*pi;
    end
    
    while (output(4) - input(4) < -pi)
        output(4) = output(4) + 2*pi;
    end

end
