function [output] = VDot(a)
    sx = [0 1; 1 0.];
    sy = [0 -1i; 1i 0.];
    sz = [1 0; 0 -1.] ;
    output = a(1)*sx + a(2)*sy + a(3)*sz;
end