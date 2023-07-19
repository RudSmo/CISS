import numpy as np
import matplotlib.pyplot as plt

import data_manager

omega = 0.1
gamma = .4
#dt_wigner = 0.1
initial_state = "DCDW"
psi_initial = "local"


wigner = data_manager.Wigner(omega, gamma, initial_state, psi_initial, load=200)


plt.pcolormesh(wigner.q_sampled, wigner.p_q_sampled, wigner.wigner, shading="gouraud")
plt.colorbar()
plt.xlabel("q")
plt.ylabel("p_q")
#plt.xlim((-4,4))
#plt.ylim((-4,4))

plt.figure()
print("mean q:", wigner.q_moment(1))
plt.plot(wigner.q_sampled, wigner.q_density())
plt.xlabel("q")
#plt.xlim((-3,3))
#plt.ylim(ymin=0)

plt.figure()
print("mean p_q:", wigner.pq_moment(1))
plt.plot(wigner.p_q_sampled, wigner.pq_density())
plt.xlabel("pq")
#plt.xlim((-3,3))
#plt.ylim(ymin=0)


plt.show()

