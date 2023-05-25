function [H] = Holstein(N,gPh,ws,xs)
    H = zeros(2*N,2*N);
    for k = 1:N
    	i = 2*k-1  ;  
	    H(i,i) = gPh(1) * sqrt(2) * xs(k);
        H(i+1,i+1) = gPh(2) * sqrt(2) * xs(k) ;
    end
  
end