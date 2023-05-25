function [output] = SOC_A_nn(N,tso,a,c,num_laps)
    Hs = zeros(2*N,2*N);
    for k = 1:N-2
	      j = 2*k-1;
        rj = UnitVector(a,c,k,num_laps,N);
        rj1 = UnitVector(a,c,k+1,num_laps,N);
        rj2 = UnitVector(a,c,k+2,num_laps,N) ;
        rjrj1 = rj-rj1;
        rjrj2 = rj-rj2;

        dj1 = rjrj1/absv(rjrj1);
        dj2 = rjrj2/absv(rjrj2);
        
        vj = cross(dj1,dj2);
         
        Vj = 1i*tso*VDot(vj);

        Hs(j,j+4) = Vj(1,1);
        Hs(j+1,j+4) = Vj(2,1);
        Hs(j,j+5) = Vj(1,2);
        Hs(j+1,j+5) = Vj(2,2) ; 
    end
    output = Hs+Hs';
end