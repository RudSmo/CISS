import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft

dat = np.genfromtxt("test_2_norm_le.dat")[:,2]
fd = fft(dat)

plt.plot(np.abs(fd))
plt.show()
