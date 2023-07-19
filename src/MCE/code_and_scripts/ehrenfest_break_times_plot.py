import numpy as np
import matplotlib.pyplot as plt

#[:, [omega, gamma, t_break]]
data = np.loadtxt("BCDW_ehrenfest_L2_break_times_num_samples=10000.dat")

plt.title("real data")
plt.scatter(data[:, 0], data[:, 1], c=np.log10(np.abs(data[:, 2])), cmap="viridis")
plt.xlabel("$\\omega$")
plt.ylabel("$\\gamma$")
cb = plt.colorbar()
cb.set_label("log10(t break)")


plt.figure()
omega_unique = np.unique(data[:, 0])
for w in omega_unique:
	which = data[:, 0] == w
	plt.plot(data[which, 1], np.abs(data[which, 2])*np.sqrt(w)*data[which, 1]**1.4, '.-', label="$\\omega=%g$"%w)
plt.xlabel("$\\gamma$")
plt.ylabel("t break*$\\sqrt{\\omega} \\gamma^{1.4}$")
plt.xlim(xmin=0)
plt.ylim(ymin=0)
plt.legend()

plt.figure()
omega_unique = np.unique(data[:, 0])
for w in omega_unique:
	which = data[:, 0] == w
	plt.plot(data[which, 1], np.abs(data[which, 2]), '.-', label="$\\omega=%g$"%w)
plt.xlabel("$\\gamma$")
plt.ylabel("t break")
plt.xlim(xmin=0)
plt.ylim(ymin=0)
plt.legend()


plt.figure()
gamma_unique = np.unique(data[:, 1])
for g in gamma_unique:
	which = data[:, 1] == g
	plt.plot(data[which, 0], np.abs(data[which, 2]), '.-', label="$\\gamma=%g$"%g)
plt.xlabel("$\\omega$")
plt.ylabel("t break")
plt.legend()




def t_break_predict(x, amp, pow_omega, pow_gamma):
	return amp*x[0]**pow_omega*x[1]**pow_omega

"""
def jac_t_break_predict(x, amp, pow_omega, pow_gamma):
	jac = np.empty((
	return amp*x[0]**pow_omega*x[1]**pow_omega

from scipy.optimize import curve_fit

p_opt, pcov = curve_fit(t_break_predict, data[:, :2].T, np.abs(data[:, 2]), p0=[1, -.5, -1.])
"""
p_opt = (1, -.5, -2.)
print(p_opt)
t_break_pred = t_break_predict(data[:, :2].T, *p_opt)

plt.figure()
plt.title("fitted modell")
plt.scatter(data[:, 0], data[:, 1], c=np.log10(t_break_pred), cmap="viridis")
plt.xlabel("$\\omega$")
plt.ylabel("$\\gamma$")
cb = plt.colorbar()
cb.set_label("log10(t break)")
plt.show()
