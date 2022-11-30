import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from matplotlib import rc
from pylab import rcParams
import random
plt.style.use('classic')

ls = 17
rc('font', family='serif', serif='Times')
rc('text', usetex=True)
rc('xtick', labelsize=ls)
rc('ytick', labelsize=ls)
rc('legend', fontsize=ls)
rc('axes', labelsize=ls)


beta = 1/0.025
w = 0.05
aj = np.tanh(0.5 * beta * w)

X = np.linspace(-15,15,5000)
W = aj/(2*np.pi) * np.exp(- aj * w * X ** 2)

Wn = random.choices(X,W,k=500)

plt.plot(X,W,color="red",linewidth=3); 
plt.hist(Wn,bins=20,density=True,color="k")
#plt.show()
plt.hist( [np.random.randn()/np.sqrt(aj*w) for i in range(9000)],bins=20,density=True,color="orange",alpha=0.3); 
plt.show()

#plt.hist()
#plt.ylabel(r"$p(X,P=0)$")
#plt.xlabel(r"$X$")
#plt.xlim(-15,15)
#plt.savefig("Comparison_initial_sampling_fake.pdf")
