import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy import constants as sc
plt.style.use("classic")

rc("font",family="serif",serif="Times")
rc("text",usetex=True)


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


fig = plt.figure(figsize=(9,12))
labels = [r"$S_0 S_0$", r"$T_1 T_1$", r"$S_1 S_0$", r"$S_0 S_1$", r"$CA$", r"$AC$", r"$DE_1$", r"$DE_2$"]


px = 0.08
py = 0.77
sx = 0.38
sy = 0.22
dx = 0.13
dy = 0.02
lw = 3
datLSC = np.load("Averaged_pops_100000_spinLSC.npy",allow_pickle=True)[()]
datPy = np.load("Averaged_pops_995000_python_dt5.npy",allow_pickle=True)[()]
datE = np.genfromtxt("exact-ortho-28.dat")

for i in range(4):
	
	for j in range(2):
		if j == 0:
			ind = i
		if j == 1:
			ind = 4+i
		pxn = px+j*(sx+dx)	
		pyn = py-i*(sy+dy)
		ax = plt.axes([pxn,pyn,sx,sy])
		ax.plot(Time(datLSC["time"],"au").get("fs"),datLSC[f"{ind}"],label="spin LSC",linewidth=lw)
		ax.plot(Time(datPy["time"],"au").get("fs"),datPy[f"2{ind}"]/64,label="spin PLDM, python",linewidth=lw)
		ax.plot(datE[:,0],datE[:,ind+1],label="MCTDH",linewidth=lw)
		ax.set_ylabel(f"Population {labels[ind]}")
			
		if i == 0 and j == 0:
			ax.legend()
		if i != 3:
			ax.set_xticklabels([])
		else:
			ax.set_xlabel(r"$t\,[\mathrm{fs}]$")
		



plt.savefig("SFOrtho_PLDM_LSC_exact_comparison.pdf")
