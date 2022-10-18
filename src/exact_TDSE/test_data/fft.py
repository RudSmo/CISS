import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft

dat = np.genfromtxt("pops_norm_le.dat")[:,2]+1j*np.genfromtxt("pops_norm_le.dat")[:,3]
fd = fft(dat)

plt.plot(np.abs(fd))
plt.show()
