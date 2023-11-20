using DelimitedFiles

tA = 50000
dat_avg = readdlm("L8_tPh_2.0_thPd_2.0_SOC0.5_SOCPh0.05_V240.0_w00.002_P0.5_SymplecticRuth4Step_RK4_VNL!_2_dt0.0002_1_currents.dat")
N = length(dat_avg[tA:end,1])
cur_avg = 1/ N * sum(dat_avg[tA:end,2:end],dims=1)
cur_std = sqrt.(1/(N-1) .* sum((cur_avg .- dat_avg[tA:end,2:end]) .^ 2,dims=1))/sqrt(N)
fil = open("currents_time_avg.dat","w")
for i in cur_avg
	write(fil,"$(i) ")
end
write(fil,"\n")
for i in cur_std
	write(fil,"$(i) ")
end
close(fil)
