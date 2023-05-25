using PyPlot
using DelimitedFiles
using FnMatch

Vs = 0:20:400
outfil = open("current_voltage.dat","w")
for v in Vs
        dat_cur = readdlm("L2_tPh_0.0_thPd_0.0_SOC0.5_SOCPh0.0_V$(v).0_w00.002_P0.5_SymplecticEulerStep_RK4_VNL!_dt0.0002_1_currents.dat")
	write(outfil,"$(v) $(dat_cur[end,2]-dat_cur[end,3])\n")
end
close(outfil)
