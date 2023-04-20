using DelimitedFiles

Vs = [48,55,65]
nS = 51
tA = 4000 
for v in Vs
	dat_avg = readdlm("V$(v)/currents_avg_n$(nS).dat")
	N = length(dat_avg[tA:end,1])
	cur_avg = 1/ N * sum(dat_avg[tA:end,2:end],dims=1)
	dat_std = readdlm("V$(v)/currents_std_n$(nS).dat")
	cur_std = sqrt.(sum(dat_std[tA:end,2:end].^2/N^2,dims=1))
	fil = open("V$(v)/currents_time_avg.dat","w")
	for i in cur_avg
		write(fil,"$(i) ")
	end
	write(fil,"\n")
	for i in cur_std
		write(fil,"$(i) ")
	end
	close(fil)
end

	
