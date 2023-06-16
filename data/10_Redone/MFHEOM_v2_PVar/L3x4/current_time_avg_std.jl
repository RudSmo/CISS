using DelimitedFiles

Vs = [50,60,75]
Ns = [11,11,11]
tA = 25000 
for (v,nS) in zip(Vs,Ns)
	for sp in ["P","M"]
		dat_avg = readdlm("V$(v)/currents_avg_n$(nS)_$(sp).dat")
		N = length(dat_avg[tA:end,1])
		cur_avg = 1/ N * sum(dat_avg[tA:end,2:end],dims=1)
		cur_std = sqrt.(1/(N-1) .* sum((dat_avg[tA:end,2:end] .- cur_avg) .^ 2,dims=1)) 
		cur_err = cur_std ./ sqrt(N)
	#	dat_std = readdlm("V$(v)/currents_std_n$(nS)_$(sp).dat")
	#	cur_std = sqrt.(sum(dat_std[tA:end,2:end].^2/N^2,dims=1))
		fil = open("V$(v)/currents_time_avg_nS$(nS)_$(sp)_std.dat","w")
		for i in cur_avg
			write(fil,"$(i) ")
		end
		write(fil,"\n")
		for i in cur_std
			write(fil,"$(i) ")
		end
		write(fil,"\n")
		for i in cur_err
			write(fil,"$(i) ")
		end
		close(fil)
	end
end

	
