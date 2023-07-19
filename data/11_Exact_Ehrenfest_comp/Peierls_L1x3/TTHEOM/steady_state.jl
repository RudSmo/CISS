using LinearAlgebra
using DelimitedFiles
using FnMatch


Vs = [22,45,50,55,60,65,70,75,85,95,105]
fil = open("TTHEOM_Iss.dat","w")
fac = 0.2*0.04
fac_c = 1e3 * 0.243394
for V in Vs
	datname = fnfilter(readdir("V$(V)"),"current*")[1]
	dat = readdlm("V$(V)/$(datname)")
	write(fil,"$(V*fac) $((dat[end,2]-dat[end,3])*fac_c)\n")
end
close(fil)
	
		
