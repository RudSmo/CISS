function [output] = UnitVector(a,c,j,num_laps,N)
    phij = 2*pi*(j-1)*num_laps/(N-1);
    cj = c*(j-1)/(N-1);
    output = [a*cos(phij),a*sin(phij),cj];
end