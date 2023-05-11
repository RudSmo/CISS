using PyPlot
using DelimitedFiles
using FnMatch

Vs = FnMatch.fnfilter(readdir("."),"V*") 
Vs = [parse(Float64,v[2:end]) for v in Vs]
Vs = sort(Vs)
outfil = open("conv_voltage.dat","w")
for v in Vs
        dat_cur = readdlm("V$(v)/"*FnMatch.fnfilter(readdir("V$(v)/"),"*convergence*")[2])
	write(outfil,"$(v) $(dat_cur[end,2]) $(dat_cur[end,3]) $(dat_cur[end,4])\n")
end
close(outfil)
