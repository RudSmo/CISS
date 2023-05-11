using PyPlot
using DelimitedFiles
using FnMatch

Vs = FnMatch.fnfilter(readdir("."),"V*") 
Vs = [parse(Float64,v[2:end]) for v in Vs]
Vs = sort(Vs)
outfil = open("current_voltage.dat","w")
for v in Vs
        dat_cur = readdlm("V$(v)/"*FnMatch.fnfilter(readdir("V$(v)/"),"*currents*")[2])
	write(outfil,"$(v) $(dat_cur[end,2]-dat_cur[end,3])\n")
end
close(outfil)
