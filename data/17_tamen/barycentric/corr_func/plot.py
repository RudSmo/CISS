import numpy as np
import matplotlib.pyplot as plt

import matplotlib.pyplot as plt;
import numpy as np;
coeff_pade=[0.-0.005j,0.-0.005j,0.-0.005j,0.-0.005j,0.-0.0050019j,0.-0.00521986j,0.-0.00765266j,0.-0.01600909j,0.-0.04639594j,0.-0.42472056j]
exp_pade=[0.01570796-0.3j,0.04712389-0.3j,0.07853982-0.3j,0.10995574-0.3j,0.14137349-0.3j,0.17313208-0.3j,0.2115157 -0.3j,0.28126243-0.3j,0.45357835-0.3j,1.33948941-0.3j]
#coeff_bary=[2.53287807e-05-0.00501972j, 7.00983951e-04-0.00652676j, 2.38371943e-03-0.01444167j, 6.33875693e-03-0.03465001j, 2.49744348e-02-0.09630535j, 4.74962281e-01-0.50494549j]
#exp_bary=[0.01572315-0.2999753j,0.04998426-0.29832741j,0.11148699-0.2890839j,0.25594994-0.26405547j,0.62627202-0.18332075j,2.0542527 +0.58499349j]
coeff_bary=np.genfromtxt("Fermi_bary_coeff.dat",dtype="complex")
exp_bary=np.genfromtxt("Fermi_bary_exp.dat",dtype="complex")
coeff_pade_e = np.genfromtxt("Fermi_coeff_2000.dat",dtype="complex")
exp_pade_e=np.genfromtxt("Fermi_exp_2000.dat",dtype="complex")[:,0]
t_au2fs = 0.024188843 
trng = np.linspace(0,10000,1000) * t_au2fs
C_pade = [np.sum([coeff*np.exp(-exp*t)*(1.+0j) for coeff,exp in zip(coeff_pade,exp_pade)]) for t in trng]
C_bary = [np.sum([coeff*np.exp(-exp*t)*(1.+0j) for coeff,exp in zip(coeff_bary,exp_bary)]) for t in trng]
C_pade_e = [np.sum([coeff*np.exp(-exp*t)*(1.+0j) for coeff,exp in zip(coeff_pade_e,exp_pade_e)]) for t in trng]

plt.plot(np.real(C_bary),label="Barycentric",color="red")
plt.plot(np.real(C_pade),label="Pade",color="blue")
plt.plot(np.real(C_pade_e),label="Pade, Np=2000",color="green")

plt.plot(np.imag(C_bary),color="red",linestyle="dashed")
plt.plot(np.imag(C_pade),color="blue",linestyle="dashed")
plt.plot(np.imag(C_pade_e),color="green",linestyle="dashed")

plt.legend()
plt.ylabel(r"$\mathrm{Re}\,  C_{L,\uparrow}^-$")
plt.xlabel(r"$\mathrm{time}\ [\mathrm{eV}^{-1}]$")
plt.ylim(-0.1,0.1)
plt.show()


