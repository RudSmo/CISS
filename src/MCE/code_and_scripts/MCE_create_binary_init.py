import numpy as np
from src.ehrenfest_init import initial_state_MCE
import os

path = "MCE_binary_start/"
os.system("mkdir -p " + path)
initial_state= "DCDW"
psi_initial  = "local"
psi_init_eigen = 1

omega = 1
gamma = 2.5
hopping= 2.5


num_sites_x = 11
num_sites_y = 1 # num_sites_x

num_configs = 1500


text_init = initial_state + "_" + psi_initial + (str(psi_init_eigen) if psi_init_eigen == "eigen" else "")
file_name = path+"MCE_%s_Lx=%lu_Ly=%lu_hopping=%g_omega=%g_gamma=%g_numconf=%i.bin"\
                  %(text_init, num_sites_x, num_sites_y, hopping, omega, gamma, num_configs)
file_raw_start = open(file_name, "wb")

hoz, el_coeff = initial_state_MCE(num_sites_x, num_sites_y, num_configs, num_configs,
                                  initial_state, psi_initial,
                                  hopping, omega, gamma, psi_init_eigen=psi_init_eigen)

print("saving data of size %i Bytes to %s"%(2*num_sites_x*num_sites_y*num_configs*16, file_name))
hoz.tofile(file_raw_start)
el_coeff.tofile(file_raw_start)
file_raw_start.close()
quit()
