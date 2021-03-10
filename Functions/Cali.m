function [output] = Cali(output,input)
    % 5G mmWave Simultaneous Localization and Mapping
    % (c) Yu Ge, 2021 (Ph.D. student at Chalmers Univerisy of Technology, emai: yuge@chalmers.se)
    % Usage: this code generates calibrated measurement
    % The measurement of the azimuth should consider its period 2pi
    % input: before calibration, measurement
    % output: after cliabration, measurement
    
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
