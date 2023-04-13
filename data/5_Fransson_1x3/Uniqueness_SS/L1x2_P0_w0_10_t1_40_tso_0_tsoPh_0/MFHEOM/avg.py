import numpy as np
import fnmatch,os
from scipy.signal import find_peaks

m = "w"
nS = 15
observables = ["currents","electrons","electrons_es","phonons"]
p1 = f"data1/"
for obs in observables:
	dat_tmp = np.genfromtxt(p1+np.array(fnmatch.filter(os.listdir(p1),f"*{obs}*"))[0])
	dat = np.zeros(np.shape(dat_tmp))
	dat[:,0] = dat_tmp[:,0]
	for j in range(2,nS+1):
	    pj = f"data{j}/"
	    dat[:,1::] += np.genfromtxt(pj+np.array(fnmatch.filter(os.listdir(pj),f"*{obs}*"))[0])[:,1::]

	dat[:,1::] = dat[:,1::]/(nS-1)
	np.savetxt(f"{obs}_avg_n{nS}.dat",dat,delimiter="\t",newline="\n")
