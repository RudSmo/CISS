function [H] = Onsite_A(N,E)
    eu = E(1);
    ed = E(2);
    
    H = zeros(2*N,2*N);
    for i = 1:2:2*N-1
        H(i,i) = eu ;
    end
    for i = 2:2:2*N
        H(i,i) = ed ;
    end

end