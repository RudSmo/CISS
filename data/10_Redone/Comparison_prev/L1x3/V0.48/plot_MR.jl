using PyPlot
using LinearAlgebra
using DelimitedFiles

plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 17
lw = 3
fac_c = 1e3 * 0.2433946
fac_c_SS = 1e6

rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)

px = 0.1
py = 0.12
sx = 0.87
sy = 0.85
dy = 0.03

fig = plt.figure(figsize=(9,5))
ax1 = plt.axes([px,py,sx,sy])
V = 0.48
colors = ["black","blue","red"]

ws = 0.002
T0 = 0.04/ws

fac_c_v2 = 1e3 * 0.2433946
fac_t_v2 = 1/(4.135668)
fac_c_v3 = 0.002 * 1e3 * 0.2433946
fac_t_v3 = 1/(ws * 4.135668)

Gamma = T0/4
annx = 0.02
anny = 0.9
ylabels = [L"$I\,[\mu\mathrm{A}]$",L"$s_1^z$"]

datMFP = readdlm("MFHEOM_v3_PVar/currents_avg_n51_P.dat")
datMFM = readdlm("MFHEOM_v3_PVar/currents_avg_n51_M.dat")
IP = datMFP[2:end,2]-datMFP[2:end,3]
IM = datMFM[2:end,2]-datMFM[2:end,3]
ID = IP .- IM
IS = IP .+ IM
SP = 100 .* ID ./ IS

tend = 11350 
dat100P = readdlm("TTHEOM/np12_r100_nV30_V0.48_P0.5/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_100.txt")
dat100M = readdlm("TTHEOM/np12_r100_nV30_V0.48_P-0.5/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_100.txt")
T100 = dat100P[2:tend,1]
I100P = dat100P[2:tend,2]-dat100P[2:tend,3]+dat100P[2:tend,4]-dat100P[2:tend,5]
I100M = dat100M[2:tend,2]-dat100M[2:tend,3]+dat100M[2:tend,4]-dat100M[2:tend,5]
I100D = I100P - I100M
I100S = I100P + I100M
SP100 = 100 .* I100D ./ I100S

tend = 7715 
dat120P = readdlm("TTHEOM/np12_r120_nV30_V0.48_P0.5/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_120.txt")
dat120M = readdlm("TTHEOM/np12_r120_nV30_V0.48_P-0.5/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_120.txt")
T120 = dat120P[2:tend,1]
I120P = dat120P[2:tend,2]-dat120P[2:tend,3]+dat120P[2:tend,4]-dat120P[2:tend,5]
I120M = dat120M[2:tend,2]-dat120M[2:tend,3]+dat120M[2:tend,4]-dat120M[2:tend,5]
I120D = I120P - I120M
I120S = I120P + I120M
SP120 = 100 .* I120D ./ I120S

tend = 4290 
dat160P = readdlm("TTHEOM/np12_r160_nV30_V0.48_P0.5/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_160.txt")
dat160M = readdlm("TTHEOM/np12_r160_nV30_V0.48_P-0.5/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_160.txt")
T160 = dat160P[2:tend,1]
I160P = dat160P[2:tend,2]-dat160P[2:tend,3]+dat160P[2:tend,4]-dat160P[2:tend,5]
I160M = dat160M[2:tend,2]-dat160M[2:tend,3]+dat160M[2:tend,4]-dat160M[2:tend,5]
I160D = I160P - I160M
I160S = I160P + I160M
SP160 = 100 .* I160D ./ I160S

tend = 3310 
dat180P = readdlm("TTHEOM/np12_r180_nV30_V0.48_P0.5/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_180.txt")
dat180M = readdlm("TTHEOM/np12_r180_nV30_V0.48_P-0.5/current_realtime_Hop+Onsite_bias_0.48_dt_0.25_rank_180.txt")
T180 = dat180P[2:tend,1]
I180P = dat180P[2:tend,2]-dat180P[2:tend,3]+dat180P[2:tend,4]-dat180P[2:tend,5]
I180M = dat180M[2:tend,2]-dat180M[2:tend,3]+dat180M[2:tend,4]-dat180M[2:tend,5]
I180D = I180P - I180M
I180S = I180P + I180M
SP180 = 100 .* I180D ./ I180S



ax1.plot(datMFP[2:end,1]*fac_t_v3,SP,color="red",label="Ehrenfest",linewidth=lw+2)
ax1.plot(T100*fac_t_v2, SP100,color= "black", label = L"$r=100$", linewidth=lw+2)
ax1.plot(T120*fac_t_v2, SP120,color= "yellow", label = L"$r=120$", linewidth=lw+1)
ax1.plot(T160*fac_t_v2, SP160,color= "turquoise", label = L"$r=160$", linewidth=lw)
ax1.plot(T180*fac_t_v2, SP180,color= "orange", label = L"$r=180$", linewidth=lw-1)

ax1.set_ylim(-3,3)
ax1.legend(loc="upper right",ncol=2)
ax1.set_ylabel(L"$M_R\,[\%]$")
ax1.set_xlabel(L"$t\,[\mathrm{fs}]$")
ax1.set_xlim(0,600)


plt.savefig("L3_comparison_MR_TTHEOM_PVar.pdf")

