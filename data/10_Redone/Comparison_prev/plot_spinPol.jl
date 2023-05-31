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
sx = 0.25
sy = 0.9
dy = 0.03
dx = 0.03

fig = plt.figure(figsize=(9,7))
ax1 = plt.axes([px,py,sx,sy])
ax2 = plt.axes([px+sx+dx,py,sx,sy])
ax3 = plt.axes([px+2*(sx+dx),py,sx,sy])
V = 55

ws = 0.002
T0 = 0.04/ws

fac_c_v2 = 1e3 * 0.2433946
fac_t_v2 = 1/(4.135668*1000)
fac_c_v3 = 0.002 * 1e3 * 0.2433946
fac_t_v3 = 1/(ws * 4.135668*1000)

Gamma = T0/4
annx = 0.02
anny = 0.9
ylabels = [L"$I\,[\mu\mathrm{A}]$",L"$s_1^z$"]


dat24_v2 = readdlm("../../9_Fransson_IV/MFHEOM/w0_2_t1_4_tso_1_tsoPh_0.1/P0.5/L2x4/V55/electrons_avg_n1515_comb.dat")
dat24_v3 = readdlm("L2x4/V55/MFHEOM_v3/P/electrons_avg_n41.dat")
dat24_B = readdlm("bundles/V55/L2x4/electrons_avg_n5_P.dat")


ax1.imshow(dat24_v2[2:end,2:2:end]-dat24_v2[2:end,3:2:end],aspect="auto")
ax2.imshow(dat24_v3[2:end,2:2:end]-dat24_v3[2:end,3:2:end],aspect="auto")
ax3.imshow(dat24_B[2:end,2:2:end]-dat24_B[2:end,3:2:end],aspect="auto")

ax1.annotate("v1",(0.8,0.07),xycoords="axes fraction",fontsize=fs,color="white")
ax2.annotate("v2",(0.8,0.07),xycoords="axes fraction",fontsize=fs,color="white")
ax3.annotate("bundles",(0.8,0.07),xycoords="axes fraction",fontsize=fs,color="white")

plt.savefig("Comparison_MFHEOM_v3_v2_bundles_spinPol_L2x4_V55.pdf")

