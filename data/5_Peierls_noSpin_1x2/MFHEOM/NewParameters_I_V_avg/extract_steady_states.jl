using LinearAlgebra
using DelimitedFiles


fil = open("I_V_n30.dat","w")
for i in 20:2:50 
	dat = readdlm("V$(i)/currents_avg_n30.dat")
	write(fil,"$(i*0.1) $(dat[end,2]-dat[end,3]) \n")
end

	
