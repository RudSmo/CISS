import numpy as np
import matplotlib.pyplot as plt

import data_manager

omega = 0.4
gamma = 1.
initial_state = "DCDW"
psi_initial = "high"

wigner = data_manager.Wigner(omega, gamma, initial_state, psi_initial)

# ref data from david
try:
	ref = data_manager.RefData(initial_state, psi_initial, 2, omega, gamma)
	ref.set_cut_time(wigner.time[-1])
except Exception as e:
	print("no reference Data found:")
	print(repr(e))
	ref = None

fig, axs = plt.subplots(3,3, figsize=(16,9))
fig.suptitle("%s %s: $\\omega=%g$, $\\gamma=%g$"%(initial_state, psi_initial, omega, gamma))
ax_qt = axs[0,0]
ax_pt = axs[0,1]
ax_qp = axs[0,2]

ax_m2 = axs[1,0]
ax_m3 = axs[1,1]
ax_m4 = axs[1,2]

ax_vari = axs[2,0]
ax_skew = axs[2,1]
ax_kurt = axs[2,2]



ax_qt.plot(wigner.time, wigner.q1, '-', label="wigner")
if ref is not None:
	ax_qt.plot(ref.time(), ref.q(), '--', label="ED")
ax_qt.set_xlabel("time t*t0")
ax_qt.set_ylabel("q")
ax_qt.legend()

ax_pt.plot(wigner.time, wigner.pq1, '-', label="wigner")
if ref is not None:
	ax_pt.plot(ref.time(), ref.pq(), '--', label="ED")
ax_pt.set_xlabel("time t*t0")
ax_pt.set_ylabel("pq")
ax_pt.legend()

ax_qp.plot(wigner.q1, wigner.pq1, '-', label="wigner")
if ref is not None:
	ax_qp.plot(ref.q(), ref.pq(), '--', label="ED")
ax_qp.set_xlabel("q")
ax_qp.set_ylabel("pq")
ax_qp.legend()


ax_m2.plot(wigner.time, wigner.q2, '-', label="q wigner")
ax_m2.plot(wigner.time, wigner.pq2, '-', label="pq wigner")
if ref is not None:
	ax_m2.plot(ref.time(), ref.q2(), '--', label="q ED")
	ax_m2.plot(ref.time(), ref.pq2(), '--', label="pq ED")
ax_m2.set_xlabel("time t*t0")
ax_m2.set_ylabel("q^2 or pq^2")
ax_m2.legend()


ax_m3.plot(wigner.time, wigner.q3, '-', label="q wigner")
ax_m3.plot(wigner.time, wigner.pq3, '-', label="pq wigner")
if ref is not None:
	ax_m3.plot(ref.time(), ref.q3(), '--', label="q ED")
	ax_m3.plot(ref.time(), ref.pq3(), '--', label="pq ED")
ax_m3.set_xlabel("time t*t0")
ax_m3.set_ylabel("q^3 or pq^3")
ax_m3.legend()


ax_m4.plot(wigner.time, wigner.q4, '-', label="q wigner")
ax_m4.plot(wigner.time, wigner.pq4, '-', label="pq wigner")
if ref is not None:
	ax_m4.plot(ref.time(), ref.q4(), '--', label="q ED")
	ax_m4.plot(ref.time(), ref.pq4(), '--', label="pq ED")
ax_m4.set_xlabel("time t*t0")
ax_m4.set_ylabel("q^4 or pq^4")
ax_m4.legend()


var_q = wigner.q2 - wigner.q1**2
var_pq = wigner.pq2 - wigner.pq1**2

if ref is not None:
	var_q_ref  =  ref.q2() -  ref.q()**2
	var_pq_ref = ref.pq2() - ref.pq()**2


ax_vari.plot(wigner.time, var_q, '-', label="q wigner")
ax_vari.plot(wigner.time, var_pq, '-', label="pq wigner")
if ref is not None:
	ax_vari.plot(ref.time(), var_q_ref, '--', label="q ED")
	ax_vari.plot(ref.time(), var_pq_ref, '--', label="pq ED")
ax_vari.set_xlabel("time t*t0")
ax_vari.set_ylabel("variance q or pq")
ax_vari.legend()


ax_skew.plot(wigner.time, ( wigner.q3 - 3* wigner.q1* wigner.q2 + 2* wigner.q1**3)/var_q**1.5, '-', label="q wigner")
ax_skew.plot(wigner.time, (wigner.pq3 - 3*wigner.pq1*wigner.pq2 + 2*wigner.pq1**3)/var_pq**1.5, '-', label="pq wigner")
if ref is not None:
	ax_skew.plot(ref.time(), ( ref.q3() - 3* ref.q()* ref.q2() + 2* ref.q()**3)/var_q_ref**1.5, '--', label="q ED")
	ax_skew.plot(ref.time(), (ref.pq3() - 3*ref.pq()*ref.pq2() + 2*ref.pq()**3)/var_pq_ref**1.5, '--', label="pq ED")
ax_skew.set_xlabel("time t*t0")
ax_skew.set_ylabel("skew q or pq")
ax_skew.legend()


ax_kurt.plot(wigner.time, ( wigner.q4 - 4* wigner.q1* wigner.q3 + 6* wigner.q1**2* wigner.q2-3* wigner.q1**4)/(var_q**2) - 3, '-', label="q wigner")
ax_kurt.plot(wigner.time, (wigner.pq4 - 4*wigner.pq1*wigner.pq3 + 6*wigner.pq1**2*wigner.pq2-3*wigner.pq1**4)/(var_pq**2) - 3, '-', label="pq wigner")
if ref is not None:
	ax_kurt.plot(ref.time(), ( ref.q4() - 4* ref.q()* ref.q3() + 6* ref.q()**2* ref.q2()-3* ref.q()**4)/(var_q_ref**2) - 3, '--', label="q ED")
	ax_kurt.plot(ref.time(), (ref.pq4() - 4*ref.pq()*ref.pq3() + 6*ref.pq()**2*ref.pq2()-3*ref.pq()**4)/(var_pq_ref**2) - 3, '--', label="pq ED")
ax_kurt.set_xlabel("time t*t0")
ax_kurt.set_ylabel("excess kurtosis q or pq")
ax_kurt.legend()

fig.tight_layout()

fig.savefig(wigner.path+"q_pq_moments.png", dpi=400)

plt.figure()
plt.title("%s %s: $\\omega=%g$, $\\gamma=%g$"%(initial_state, psi_initial, omega, gamma))
plt.xlabel("time")
plt.ylabel(r"$\langle n \rangle$")
plt.tight_layout()
try:
	plt.plot(wigner.time, wigner.nel())
	plt.savefig(wigner.path + "electron_density.png")
except FileNotFoundError:
	pass

print(wigner.q1[0], wigner.pq1[0])
print(var_q[0], var_pq[0])
plt.show()
