import numpy as np
import fnmatch,os

m = "w"
observables = ["currents","electrons","electrons_es","phonons","convergence_info"]
Ns = 11 
for obs in observables:
	p1 = f"data0/"
	dat_tmp = np.genfromtxt(p1+np.array(fnmatch.filter(os.listdir(p1),f"*{obs}_avg.dat"))[0])
	dat = np.zeros(np.shape(dat_tmp))
	dat[:,0] = dat_tmp[:,0]
	for j in range(1,Ns):
		pj = f"data{j}/"
		dat[:,1::] += np.genfromtxt(pj+np.array(fnmatch.filter(os.listdir(pj),f"*{obs}_avg.dat"))[0])[:,1::]

	dat[:,1::] = dat[:,1::]/(Ns)
	np.savetxt(f"{obs}_avg_n{Ns*3000}.dat",dat,delimiter="\t",newline="\n")
