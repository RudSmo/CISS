using LinearAlgebra
using DelimitedFiles

fac_c =1e3 * 0.2433946 
V = [0,10,11,12,13,14,15,16,17,18,25]
w0 = 0.05
t1 = 0.02 

fil = open("I_V.dat","w")

for v in V
	fname = "V$(v)/"
	
	ve = v * 0.04
	curTT = readdlm("$(fname)current_realtime_Hop+Onsite_bias_$(ve)_dt_0.25_rank_100.txt")
	write(fil,"$(v) $((curTT[end,2]-curTT[end,3])*fac_c) \n")
end
close(fil)

