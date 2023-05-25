
function [thing] = FranssonHamiltonian(N,OnsiteEnergies,Hopping,SOC,HoppingPh,SOCPh,a_hel,c_hel,num_laps,ws,X)
	Ho = Onsite_A(N,OnsiteEnergies);
	Hh = Hopping_A_nn(N,Hopping);
	HSo = SOC_A_nn(N,SOC,a_hel,c_hel,num_laps);
  Heph = Peierls(N,HoppingPh,ws,X);
  Hsoph = Peierls_SOC(N,SOCPh,a_hel,c_hel,num_laps,ws,X);
      
  thing = Ho + Hh + HSo + Heph + Hsoph;
  
end

