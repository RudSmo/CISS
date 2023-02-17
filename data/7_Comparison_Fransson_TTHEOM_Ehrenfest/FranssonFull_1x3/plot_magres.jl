using LinearAlgebra
using PyPlot
using DelimitedFiles

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

w0 = 0.05
V = 0.48
px = 0.1
py = 0.13
sx = 0.875
sy = 0.85
dy = 0.05

fig = plt.figure(figsize=(10,5))
ax1 = plt.axes([px,py,sx,sy])


curMF1 = readdlm("MFHEOM/P0.5/currents_avg_n30.dat")
curMF2 = readdlm("MFHEOM/P-0.5/currents_avg_n30.dat")

curNo1 = readdlm("noPh/L3_tPh0.0_SOC0.001_SOCPh0.0_V0.48_w00.002_P0.5_X0_0.0_dt0.05_1_currents.dat") 
curNo2 = readdlm("noPh/L3_tPh0.0_SOC0.001_SOCPh0.0_V0.48_w00.002_P-0.5_X0_0.0_dt0.05_1_currents.dat")

t_end = 1000

P = 0.5
TT_Cur_Files1 = ["TTHEOM/P$(P)/1x3_sym_np10_r100_nV30_V0.48_P$(P)_corrected/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_100.txt",
                "TTHEOM/P$(P)/1x3_sym_np10_r120_nV30_V0.48_P$(P)_corrected/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_120.txt",
                "TTHEOM/P$(P)/1x3_sym_np10_r140_nV30_V0.48_P$(P)_corrected/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_140.txt"]
TT_Cur_Files2 = ["TTHEOM/P-$(P)/1x3_sym_np10_r100_nV30_V0.48_P-$(P)_corrected/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_100.txt",
                "TTHEOM/P-$(P)/1x3_sym_np10_r120_nV30_V0.48_P-$(P)_corrected/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_120.txt",
                "TTHEOM/P-$(P)/1x3_sym_np10_r140_nV30_V0.48_P-$(P)_corrected/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_140.txt"]

TT_labels = [L"$r = 100$",
             L"$r = 120$",
             L"$r = 140$"]


I1 = (curMF1[2:end,2]-curMF1[2:end,3])*fac_c*2
I2 = (curMF2[2:end,2]-curMF2[2:end,3])*fac_c*2
M_MF = (I1 .- I2) ./ (I1 .+ I2)

I1no = (curNo1[2:end,2]-curNo1[2:end,3])*fac_c*2
I2no = (curNo2[2:end,2]-curNo2[2:end,3])*fac_c*2
M_no = (I1no .- I2no) ./ (I1no .+ I2no)

ax1.plot(curMF1[2:end,1]*fac_t, M_MF*100, color="k",label = "MF-HEOM", linewidth=lw+2)
ax1.plot(curNo1[2:end,1]*fac_t, M_no*100, color="orange",label = L"exact $t_1=\lambda_1=\, 0.0\mathrm{meV}$", linewidth=lw, linestyle="dashed")

colors = ["green", "blue", "red"]
ends = [12000,12000,8500]
for (cur_fil1,cur_fil2,label,col,a) in zip(TT_Cur_Files1,TT_Cur_Files2,TT_labels,colors,ends)
	cur_dat1 = readdlm(cur_fil1)[2:a,:]
	cur_dat2 = readdlm(cur_fil2)[2:a,:]
	ti_cur1 = cur_dat1[:,1]
	ti_cur2 = cur_dat2[:,1]
#	cur_dat1 = cur_dat1[sortperm(ti_cur1),:]
#	cur_dat2 = cur_dat2[sortperm(ti_cur2),:]
	
	I1t = (cur_dat1[:,2]-cur_dat1[:,3]+cur_dat1[:,4]-cur_dat1[:,5])*fac_c	
	I2t = (cur_dat2[:,2]-cur_dat2[:,3]+cur_dat2[:,4]-cur_dat2[:,5])*fac_c
	Mt = (I1t .- I2t) ./ (I1t .+ I2t)
	ax1.plot(ti_cur1*fac_t,Mt*100,color=col,label=label,linewidth=lw)
end
ax1.legend(loc="lower left",numpoints=1)
ax1.set_ylabel(L"$M = \frac{I(P=0.5) -I(P=-0.5)}{I(P=0.5) +I(P=-0.5)}\,[\%]$")
ax1.set_xlabel(L"$t\,[\mathrm{fs}]$")
ax1.set_xlim(0,t_end)
ax1.set_ylim(-10,5)



plt.savefig("Comparison_TTHEOM_MFHEOM_Fransson_1x3_w0_$(w0)_V$(V)_Magnetoresistance.pdf")

