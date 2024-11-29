import matplotlib.pyplot as plt
import numpy as np

from matplotlib.collections import LineCollection
from matplotlib.colors import BoundaryNorm, ListedColormap

plt.style.use("classic")
plt.rc("text", usetex=True)
plt.rc("font",family="serif",serif="Times")
fs = 21
lw = 3
fac_c = 0.002*0.2*1e3 * 0.2433946
fac_c_SS = 1e6
fac_t = 1/(0.002*0.2*4.135668*1e3)
fac_e = 0.002*0.2

plt.rc("axes",labelsize=fs)
plt.rc("xtick",labelsize=fs-2)
plt.rc("ytick",labelsize=fs-2)

px = 0.14
py = 0.821
sx1 = 0.85
sx2 = 0.68
sy = 0.17
dy = 0.021
dx = 0.01

fig = plt.figure(figsize=(8,14))
ax1 = plt.axes([px,py,sx1,sy])
ax2 = plt.axes([px,py-sy-dy,sx2,sy])
ax3 = plt.axes([px,py-2*(sy+dy),sx1,sy])
ax4 = plt.axes([px,py-3*(sy+dy),sx2,sy])
ax5 = plt.axes([px,py-4*(sy+dy),sx2,sy])

sp = "P"
dat_sp = np.genfromtxt("spec_P_L1x5.dat")
dat_el = np.genfromtxt("L1x5/L5_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.25_V1200.0_w00.0004_P0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_electrons_es.dat")
N = len(dat_sp[0,1::])
x = dat_sp[0::10,0]*fac_t
mn = dat_el[:,1::].min()
mx = dat_el[:,1::].max()
norm = plt.Normalize(mn, mx)
for i in range(N):
	y = dat_sp[0::10,i+1]*fac_e
	c = dat_el[0::10,i+1]
	points = np.array([x, y]).T.reshape(-1, 1, 2)
	segments = np.concatenate([points[:-1], points[1:]], axis=1)
	lc = LineCollection(segments, cmap='viridis', norm=norm)
	# Set the values used for colormapping
	lc.set_array(c)
	lc.set_linewidth(2)
	line = ax1.add_collection(lc)


cbar = fig.colorbar(line, ax=ax1)
cbar.set_label("adiabatic state pop",labelpad=10)
cbar.ax.tick_params(labelsize=fs-2)
ax1.set_xlim(segments[:,:,0].min(), segments[:,:,0].max())
ax1.set_ylim(-0.5,0.0)
ax1.set_ylabel("energy "+r"$E\,[\mathrm{eV}]$")
ax1.annotate(r"$A=+\frac{1}{2}$",(0.05,0.9),xycoords="axes fraction",fontsize=fs)


sufP="_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.25_V1200.0_w00.0004_P0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_phonons.dat"
sufM="_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.25_V1200.0_w00.0004_P-0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_phonons.dat"
sufIP="_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.25_V1200.0_w00.0004_P0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_currents.dat"
sufIM="_tPh_10.0_thPd_10.0_SOC2.5_SOCPh0.25_V1200.0_w00.0004_P-0.5_SymplecticRuth4Step_RK4_VNL!_dt5.0e-5_1_currents.dat"
pre = ["L1x5/L5"]
col = ["k"]
for p,c in zip(pre,col):
	datXP = np.genfromtxt(f"{p}{sufP}")[1::,:]
	datXM = np.genfromtxt(f"{p}{sufM}")[1::,:]
	lXP = len(datXP[:,0])
	lXM = len(datXM[:,0])
	lf = min(lXP,lXM)
	datXP = datXP[0:lf,:]
	datXM = datXM[0:lf,:]
	tim = datXP[:,0]*fac_t
	XP = datXP[:,1]+datXP[:,3]+datXP[:,5]+datXP[:,7]
	XM = datXM[:,1]+datXM[:,3]+datXM[:,5]+datXM[:,7]
	dX = XP-XM
	ax2.plot(tim,dX,color=c,linewidth=lw-1)


	datP = np.genfromtxt(f"{p}{sufIP}")[1::,:]
	datM = np.genfromtxt(f"{p}{sufIM}")[1::,:]
	lP = len(datP[:,0])
	lM = len(datM[:,0])
	lD = min(lP,lM)
	ti = datP[0:lD,0]*fac_t
	IP = (datP[0:lD,1]-datP[0:lD,2])*fac_c/2
	IM = (datM[0:lD,1]-datM[0:lD,2])*fac_c/2
	ax4.plot(ti,IP,color="red",linewidth=lw-1)
	ax4.plot(ti,IM,color="blue",linewidth=lw-1)

	Is = IP+IM
	Id = IP-IM
	SP = 100 * Id / Is
	ax5.plot(ti,SP,color=c,linewidth=lw-1)



ax2.set_xlim(0,50)
ax2.set_ylim(-10,10)
ax2.set_ylabel(r"$\overline{X}_+ - \overline{X}_-$")


sp = "P"
dat_sp = np.genfromtxt("spec_P_L1x5.dat")
dat_el = np.genfromtxt("gammaL_P_L1x5.dat")
N = len(dat_sp[0,1::])
x = dat_sp[0::10,0]*fac_t
mn = dat_el[:,1::].min()*0.4
mx = dat_el[:,1::].max()*0.2
norm = plt.Normalize(mn, mx)
for i in range(N):
        y = dat_sp[0::10,i+1]*fac_e
        c = dat_el[0::10,i+1]*0.4
        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        lc = LineCollection(segments, cmap='viridis', norm=norm)
        # Set the values used for colormapping
        lc.set_array(c)
        lc.set_linewidth(2)
        line = ax3.add_collection(lc)


cbar = fig.colorbar(line, ax=ax3)
cbar.set_label(r"$\Gamma_L(t) \, [\mathrm{meV}]$",labelpad=10)
cbar.ax.tick_params(labelsize=fs-2)
ax3.set_xlim(segments[:,:,0].min(), segments[:,:,0].max())
ax3.set_ylim(-0.5,0.0)
ax3.set_ylabel("energy "+r"$E\,[\mathrm{eV}]$")
ax4.set_ylabel(r"$I(t)\, [\mu\mathrm{A}]$")
ax5.set_ylabel(r"$S(t)\, [\%]$")
ax3.set_xticklabels([])

ax1.set_xticklabels([])
ax2.set_xticklabels([])
ax3.set_xlim(0,50)
ax4.set_xticklabels([])
ax4.set_xlim(0,50)
ax5.set_xlim(0,50)
ax5.set_xlabel("time "+r"$t\, [\mathrm{ps}]$")
plt.savefig(f"L1x5_elec_es_spec_Xdif_gamma_td_cur_{sp}_singleTraj.pdf")
