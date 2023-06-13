using LinearAlgebra
using PyPlot
using DelimitedFiles
using FnMatch


plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
lw = 3
fac_c = 1e3 * 0.2433946
fac_c_SS = 1e6
fac_t = 1/(4.135668)

rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)

px = 0.09
py = 0.12
sx = 0.89
sy = 0.85
dy = 0.05
dx = 0.01

fig = plt.figure(figsize=(9,7))
ax = plt.axes([px,py,sx,sy])
P = 0.5

V = [60]
colors = ["k","red","blue"]
size = [3,4,8,12,16]
ms = 12
for (v,vn,col) in zip(V,Vn,colors)

	dats = ["L2x5/V$(v)/currents_time_avg_nS101",
		"L3x3/V$(v)/currents_time_avg_nS201",
		"L2x4/V$(v)/currents_time_avg_nS101",
		"L2x4/V$(v)/currents_time_avg_nS101",
		"L2x3/V$(v)/currents_time_avg_nS1417", 
		"L1x3/V$(v)/currents_time_avg_nS301"]


	dats_n = ["L4x4/V$(v)/noPh/L16_tPh_0.0_thPd_0.0_SOC0.001_SOCPh0.0_V$(vn)_w00.002_P0.5_noEps_dt0.05_1_currents.dat",
		  "L3x4/V$(v)/noPh/L12_tPh_0.0_thPd_0.0_SOC0.001_SOCPh0.0_V$(vn)_w00.002_P0.5_noEps_dt0.05_1_currents.dat",
	   	  "L2x4/V$(v)/noPh/L8_tPh_0.0_thPd_0.0_SOC0.001_SOCPh0.0_V$(vn)_w00.002_P0.5_noEps_dt0.05_1_currents.dat",
	   	  "L2x2/V$(v)/noPh/L4_tPh_0.0_thPd_0.0_SOC0.001_SOCPh0.0_V$(vn)_w00.002_P0.5_noEps_dt0.05_1_currents.dat",
	   	  "L1x3/V$(v)/noPh/L3_tPh_0.0_thPd_0.0_SOC0.001_SOCPh0.0_V$(vn)_w00.002_P0.5_noEps_dt0.05_1_currents.dat"]

	IT = []
	sIT = []

	ITn = []
	for dat in dats[end:-1:1]
		data = readdlm(dat*"_std.dat")
		push!(IT,(data[1,1]-data[1,2])*fac_c)
		serr = sqrt((data[2,1]*fac_c)^2 + (data[2,2]*fac_c)^2)
		push!(sIT,serr)
	end

	for dat in dats_n[end:-1:1]
		data = readdlm(dat)
		push!(ITn,(data[end,2]-data[end,3])*fac_c)
	end


	if v == V[1]
	
		ax.errorbar(size,IT,yerr=sIT,marker="o",linewidth=2,color=col,label=L"$t_1=4\,\mathrm{meV}, \ \lambda_1 = 0.1\,\mathrm{meV}$",markersize=ms,linestyle="")
		ax.plot(size,ITn,marker="d",linewidth=0,color="red",label=L"$t_1=\lambda_1 = 0\,\mathrm{meV}$",markersize=ms)
	else
		ax.errorbar(size,IT,yerr=sIT,marker="o",linewidth=1,color=col,markersize=ms)
		ax.plot(size,ITn,marker="d",linewidth=0,color="red",markersize=ms)
	end

end

ax.set_xlabel("system size")
ax.set_ylabel(L"$I\,[\mu\mathrm{A}]$")
ax.set_xlim(1,17)
ax.legend(loc="lower left")
plt.savefig("Currents_sizeVar_P$(P)_VVar.png",dpi=200)
