import numpy as np
import matplotlib.pyplot as plt

num_sites = 11
hopping = 1
gamma = 1
num_samples = 100
x = np.linspace(-3, 4, num_samples)

H = np.zeros((num_samples, num_sites, num_sites))
diag_ind = np.arange(num_sites)

H[:, diag_ind[:-1], diag_ind[1:]] = -hopping
H[:, diag_ind[1:], diag_ind[:-1]] = -hopping
H[:, diag_ind, diag_ind] = (x**2).reshape(num_samples, 1)
H[:, num_sites//2, num_sites//2] -= gamma*x

# find eigenvalues and eigenstates
eigval, eigvecs = np.linalg.eigh(H)
eigval = eigval.reshape(num_samples, num_sites)
order = np.argsort(eigval, axis=-1)

plt.plot(x, eigval)
plt.show()

