function [H] = Peierls(N,tph,ws,xs)
    H = zeros(2*N,2*N);
    for k = 1:N-1
	i = 2*k-1;
        H(i,i+2) = tph(1) * xs(k)/sqrt(2) ;
	H(i+1,i+3) = tph(2) * xs(k)/sqrt(2);
        H(i+2,i) = tph(1) * xs(k)/sqrt(2);
	H(i+3,i+1) = tph(2) * xs(k)/sqrt(2);
    end
   
end