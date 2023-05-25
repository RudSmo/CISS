function [Hh] = Hopping_A_nn(N,t)
    Hh = zeros(2*N,2*N);
    for i = 1:2*N-2
        Hh(i,i+2) = t;
        Hh(i+2,i) = t;
    end
end