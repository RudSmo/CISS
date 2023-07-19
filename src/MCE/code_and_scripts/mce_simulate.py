import numpy as np
import data_manager

from sys import argv

if len(argv) not in [13, 14, 15, 16]:
	print("invalid number of parameters, expected: t0 omega gamma Lx Ly num_configs "\
	      "dt t_max initial_state[bare, dressed] psi_initial[local, eigen, CDW] "\
	      "calc_BO_electron_density[True, False] BO_nel_save_dist psi_state[2 ints y, x (local) or 1(eigen) or 0(CDW) inst]"\
	      " [name_of_subfolder, default=0])")
	print("expected 12, got", len(argv)-1)
	quit(1)

version = 1
hopping = float(argv[1])
omega = float(argv[2])
gamma = float(argv[3])
num_sites_x = int(argv[4])
num_sites_y = int(argv[5])

num_configs = int(argv[6])
max_num_configs = num_configs

dt = float(argv[7])
t_max = float(argv[8])

initial_state = argv[9]
psi_initial = argv[10]

if psi_initial == "local":
	psi_state = (int(argv[13]), int(argv[14]))
elif psi_initial == "eigen":
	psi_state = int(argv[13])
else:
	psi_state = None


calc_BO_electron_density = argv[11].lower() == "true"
if not calc_BO_electron_density and argv[11].lower() != "false":
	print(f"argv[11]=\"{argv[11]}\" must be true or false (case insensitive)")
	quit(1)
BO_nel_save_dist = int(argv[12])
extra_subfolder=None
if    len(argv) == 16 and psi_initial == "local"\
   or len(argv) == 15 and psi_initial == "eigen"\
   or len(argv) == 14 and psi_initial == "CDW":
	extra_subfolder = str(argv[-1])

load_sampled_state = False
force_new_simulation = True


sim = data_manager.Simulation_MCE(version, hopping, omega, gamma,
                                  num_sites_x, num_sites_y,
                                  t_max, dt,
                                  num_configs,
                                  initial_state, psi_initial, psi_state,
                                  max_num_configs, printProgress=True,
                                  force_new_simulation=force_new_simulation,
                                  load_sampled_state=load_sampled_state,
                                  store_results=True,
                                  calc_BO_electron_density=calc_BO_electron_density,
                                  BO_nel_save_dist=BO_nel_save_dist,
                                  extra_subfolder=extra_subfolder)

