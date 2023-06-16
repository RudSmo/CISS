using LinearAlgebra
using PyPlot
using DelimitedFiles
using FnMatch


plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 15
lw = 3
fac_c = 0.002* 1e3 * 0.2433946
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

fig = plt.figure(figsize=(9,4.5))
ax = plt.axes([px,py,sx,sy])
P = 0.5
# OLD DATA
V = [48,55,65]
Vn = ["0.38400000000000006"]#,"0.44","0.52"]
colors = ["orange","red","blue"]
size = [3,4,8,12,16]
ms = 7
ppath = "../../9_Fransson_IV/MFHEOM/w0_2_t1_4_tso_1_tsoPh_0.1/P0.5/"
fac_c = fac_c/0.002
#for (v,vn,col) in zip(V,Vn,colors)
#
#        dats = ["$(ppath)L4x4/V$(v)/currents_time_avg_nS101",
#                "$(ppath)L3x4/V$(v)/currents_time_avg_nS201",
#                "$(ppath)L2x4/V$(v)/currents_time_avg_nS101",
#                "$(ppath)L2x2/V$(v)/currents_time_avg_nS1417",
#                "$(ppath)L1x3/V$(v)/currents_time_avg_nS301"]
#
#
#        dats_n = ["$(ppath)L4x4/V$(v)/noPh/L16_tPh_0.0_thPd_0.0_SOC0.001_SOCPh0.0_V$(vn)_w00.002_P0.5_noEps_dt0.05_1_currents.dat",
#                  "$(ppath)L3x4/V$(v)/noPh/L12_tPh_0.0_thPd_0.0_SOC0.001_SOCPh0.0_V$(vn)_w00.002_P0.5_noEps_dt0.05_1_currents.dat",
#                  "$(ppath)L2x4/V$(v)/noPh/L8_tPh_0.0_thPd_0.0_SOC0.001_SOCPh0.0_V$(vn)_w00.002_P0.5_noEps_dt0.05_1_currents.dat",
#                  "$(ppath)L2x2/V$(v)/noPh/L4_tPh_0.0_thPd_0.0_SOC0.001_SOCPh0.0_V$(vn)_w00.002_P0.5_noEps_dt0.05_1_currents.dat",
#                  "$(ppath)L1x3/V$(v)/noPh/L3_tPh_0.0_thPd_0.0_SOC0.001_SOCPh0.0_V$(vn)_w00.002_P0.5_noEps_dt0.05_1_currents.dat"]
#
#        IT = []
#        sIT = []
#
#        ITn = []
#        for dat in dats[end:-1:1]
#                data = readdlm(dat*"_std.dat")
#                push!(IT,(data[1,1]-data[1,2])*fac_c)
#                serr = sqrt((data[2,1]*fac_c)^2 + (data[2,2]*fac_c)^2)
#                push!(sIT,serr)
#        end
#
#        for dat in dats_n[end:-1:1]
#                data = readdlm(dat)
#                push!(ITn,(data[end,2]-data[end,3])*fac_c)
#        end
#
#
#        if v == V[1]
#
#                ax.errorbar(size,IT .* sqrt.(size),yerr=sIT.* sqrt.(size),marker="o",linewidth=2,color=col,label=L"$0.348\,\mathrm{V}$"*" (prev. code)",markersize=ms)
#                ax.plot(size,ITn .* sqrt.(size)/v,marker="d",linewidth=2,color="orange",label=L"$0.348\,\mathrm{V},\ t_1=\lambda_1 = 0\,\mathrm{meV}$",markersize=ms)
#        else
#                ax.errorbar(size,IT,yerr=sIT,marker="o",linewidth=1,color=col,markersize=ms)
#                ax.plot(size,ITn,marker="d",linewidth=0,color="orange",markersize=ms)
#        end
#
#end



# NEW DATA
fac_c = 0.002* 1e3 * 0.2433946
V = [50,60,75]
Vn = [L"$0.4\,\mathrm{V}$",L"$0.48\,\mathrm{V}$",L"$0.6\,\mathrm{V}$"]
Ns =[[11,11,21,21,72,10,71],
     [11,11,19,21,71,10,61],
     [11,11,20,11,72,10,61]]
colors = ["k","red","blue"]
size = [3,4,6,8,9,10,12]
ms = 8



for (i,v,vn,col) in zip(1:length(V),V,Vn,colors)

	dats = ["L3x4/V$(v)/currents_time_avg_nS$(Ns[i][1])_P",
		"L2x5/V$(v)/currents_time_avg_nS$(Ns[i][2])_P",
		"L3x3/V$(v)/currents_time_avg_nS$(Ns[i][3])_P",
		"L2x4/V$(v)/currents_time_avg_nS$(Ns[i][4])_P",
		"L2x3/V$(v)/currents_time_avg_nS$(Ns[i][5])_P",
		"L2x2/V$(v)/currents_time_avg_nS$(Ns[i][6])_P", 
		"L1x3/V$(v)/currents_time_avg_nS$(Ns[i][7])_P"]


#	dats_n = ["L4x4/V$(v)/noPh/L16_tPh_0.0_thPd_0.0_SOC0.001_SOCPh0.0_V$(vn)_w00.002_P0.5_noEps_dt0.05_1_currents.dat",
#		  "L3x4/V$(v)/noPh/L12_tPh_0.0_thPd_0.0_SOC0.001_SOCPh0.0_V$(vn)_w00.002_P0.5_noEps_dt0.05_1_currents.dat",
#	   	  "L2x4/V$(v)/noPh/L8_tPh_0.0_thPd_0.0_SOC0.001_SOCPh0.0_V$(vn)_w00.002_P0.5_noEps_dt0.05_1_currents.dat",
#	   	  "L2x2/V$(v)/noPh/L4_tPh_0.0_thPd_0.0_SOC0.001_SOCPh0.0_V$(vn)_w00.002_P0.5_noEps_dt0.05_1_currents.dat",
#	   	  "L1x3/V$(v)/noPh/L3_tPh_0.0_thPd_0.0_SOC0.001_SOCPh0.0_V$(vn)_w00.002_P0.5_noEps_dt0.05_1_currents.dat"]

	IT = []
	sIT = []

	ITn = []
	for dat in dats[end:-1:1]
		data = readdlm(dat*"_std.dat")
		push!(IT,(data[1,1]-data[1,2])*fac_c)
		serr = sqrt((data[2,1]*fac_c)^2 + (data[2,2]*fac_c)^2)
		push!(sIT,serr)
	end

#	for dat in dats_n[end:-1:1]
#		data = readdlm(dat)
#		push!(ITn,(data[end,2]-data[end,3])*fac_c)
#	end


	ax.errorbar(size,IT.* sqrt.(size),yerr=sIT .* sqrt.(size),marker="o",linewidth=2,color=col,label=vn,markersize=ms)

end
IT = []
for i in size
	dat = readdlm("Gamma_NVar/L$(i)_tPh_0.0_thPd_0.0_SOC0.5_SOCPh0.0_V240.0_w00.002_P0.5_SymplecticRuth4Step_RK4_VNL!_dt0.0002_1_currents.dat")
	push!(IT,(dat[end,2]-dat[end,3])*fac_c*sqrt(i))
end

ax.plot(size,IT,marker="o",linewidth=2,markersize=ms,label=L"$0.48\,\mathrm{V},\ \Gamma=\Gamma_0/\sqrt{N},\ t_1=\lambda_1=0$",color="green")


ax.set_xlabel("system size "*L"$N$")
ax.set_ylabel(L"$\sqrt{N}I\,[\mu\mathrm{A}]$")
ax.set_xlim(2,13)
ax.legend(loc="best",bbox_to_anchor=(0.55, 0.4, 0.5, 0.5))
#plt.savefig("Currents_sizeVar_P$(P)_VVar.png",dpi=200)
plt.savefig("Currents_sizeVar_P$(P)_VVar_sqrtN.pdf")
