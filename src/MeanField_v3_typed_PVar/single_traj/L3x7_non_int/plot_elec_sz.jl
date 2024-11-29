using PyPlot
using DelimitedFiles


plt.style.use("classic")
rc("text", usetex=true)
rc("font",family="serif",serif="Times")
fs = 19
lw = 3
fac_c = 0.002*0.2*1e3 * 0.2433946
fac_c_SS = 1e6
fac_t = 1/(0.002*0.2*4.135668*1e3)

rc("axes",labelsize=fs)
rc("xtick",labelsize=fs-2)
rc("ytick",labelsize=fs-2)

px = 0.12
py = 0.12
sx = 0.7
sy = 0.85
dy = 0.05
dx = 0.01

fig = plt.figure(figsize=(8,5))
ax = plt.axes([px,py,sx,sy])
cax = plt.axes([px+(sx+dx),py,0.03,sy])
cp = plt.cm.plasma
cmap = plt.cm.ScalarMappable(cmap=cp)

sp = "P"
dat = "electrons_avg_n1001_$(sp).dat"

N = 21
datE = readdlm("L21_tPh_0.0_thPd_0.0_SOC2.5_SOCPh0.0_V1200.0_w00.0004_P0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_electrons.dat")
SP = datE[2:end,2:2:end] .- datE[2:end,3:2:end]
tim = datE[end,1]*fac_t
im = ax.imshow(SP[end:-1:1,:],aspect="auto",vmin=0.01,vmax=0.2,extent=[1,21,0,tim],cmap="plasma")
cbar = fig.colorbar(im, ax=ax, cax = cax, aspect=90,pad=.001)
cbar.set_label("electronic spin polarization "*L"$s_{j}^z\,[\frac{\hbar}{2}]$")
cbar.ax.tick_params(labelsize=fs-2)
ax.annotate(L"$t_1=\lambda_1=0$",(0.03,0.9),xycoords="axes fraction",fontsize=fs,color="white")
ax.annotate(L"$A=+\frac{1}{2}$",(0.8,0.9),xycoords="axes fraction",fontsize=fs,color="white")
ax.set_xlabel("site "*L"$j$")
ax.set_ylabel("time "*L"$t\,[\mathrm{ps}]$",labelpad=7)
plt.savefig("L3x7_V60_t1_0_l1_0_elec_sz_$(sp).pdf")
