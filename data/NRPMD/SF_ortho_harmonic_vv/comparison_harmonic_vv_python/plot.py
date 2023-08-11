import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import cm
from scipy import constants as sc


class Quantity(object):
	def __init__(self, value, unit=''):
		self.value = value
		self.unit = unit
	def __str__(self):
		return "%g %s" % (self.value, self.unit)
	def set(self, value, unit=None):
		"""Set value in specified unit (default: same as before)."""
		self.value = value
		if unit is not None:
			self.unit = unit
	def get(self, unit=None):
		"""Return value in specified unit (default: same as before)."""
		if unit is None:
			return self.value
		else:
			value = self.value * self.U[self.unit]
			return value / self.U[unit]
	def change(self, unit):
		"""Return Quantity in new unit"""
		newval = self.get(unit)
		return Quantity(newval, unit)
class Time(Quantity):
	U = {'s':1, 'fs':sc.femto, 'ps':sc.pico, 'ns':sc.nano,'mus':sc.micro, 'au':sc.value("atomic unit of time")}
	def __init__(self, value, unit='s'):
		Quantity.__init__(self, value, unit)

plt.style.use("classic")

rc("font",family="serif",serif="Times")
rc("text",usetex=True)



dat1 = np.load("Averaged_pops_1000000_harmonic_vv_dt5.npy",allow_pickle=True)[()]
dat2 = np.load("Averaged_Pops_1000000_harmonic_vv_dt5_initialState_02.npy",allow_pickle=True)[()]
dat3 = np.load("Averaged_pops_995000_python_dt5.npy",allow_pickle=True)[()]
dat4 = np.load("Averaged_pops_285000_python_dt5_initialState_02.npy",allow_pickle=True)[()]
datE = np.genfromtxt("exact-ortho-28.dat")

dat = [dat1,dat2,dat3,dat4]

lab = ["harmonic, "+r"$n_\mathrm{traj}=10^6$","harmonic, "+"$n_\mathrm{traj}=10^6$, -initial\_state 2", "python, "+"$n_\mathrm{traj}=5\cdot 10^5$", "python, "+"$n_\mathrm{traj}=10^5$, -initial\_state 2"]

llab = [r"$S_0 S_0$", r"$T_1 T_1$", r"$S_1 S_0$", r"$S_0 S_1$", r"$CA$", r"$AC$", r"$DE_1$", r"$DE_2$"]

fig = plt.figure(figsize=(6.5,9))

px = 0.1
py = 0.765
sx = 0.8
sy = 0.215
dy = 0.02
ax0 = plt.axes([px,py,sx,sy])
ax1 = plt.axes([px,py-sy-dy,sx,sy])
ax2 = plt.axes([px,py-2*(sy+dy),sx,sy])
ax3 = plt.axes([px,py-3*(sy+dy),sx,sy])
ax = [ax0,ax1,ax2,ax3]
colors = cm.jet(np.arange(0.0,1.0,1/8))
for i,j in enumerate(dat):
	if i != 3:
		ax[i].tick_params(labelbottom=False)
	ax[i].annotate(lab[i],(0.1,0.9),xycoords="axes fraction")
	ax[i].set_ylim(-0.05,1)
	for k,col in enumerate(colors):
		li = "2"+str(k)
		time_end = j["time"][-1]
		ax[i].plot(Time(j["time"],"au").get("fs"),np.abs(j[li])/64,label=llab[k],color=col,linewidth=2)
		ax[i].set_ylabel(r"Population")

colors=cm.jet(np.arange(0.0,1.0,1/8)-0.05)
for i,c in zip(range(1,len(datE[0,:])),colors):
	ax[-1].plot(datE[:,0],datE[:,i],linestyle="--",dashes=(5, 1),color=c)



ax[0].legend(title=r"State",handlelength=1.2,fontsize=12,bbox_to_anchor=(1.13,1.1),labelspacing=0.3)
ax[-1].set_xlabel(r"Time (fs)")

plt.savefig("SFOrtho_Comparison_spin_PLDM_python_harmonic_is2.pdf")
