using LinearAlgebra
using DelimitedFiles

fac_c =1e3 * 0.2433946 
V = [0,10,11,12,13,14,15,16,17,18]
w0 = 0.05
t1 = 0.02 

fil = open("I_V.dat","w")

for v in V
	fname = "V$(v)/"
	cur = readdlm("$(fname)currents_avg_n30.dat")
	write(fil,"$(v) $((cur[end,2]-cur[end,3])*fac_c) \n")
end
close(fil)

