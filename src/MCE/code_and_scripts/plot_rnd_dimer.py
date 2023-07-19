import numpy as np
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerLine2D, HandlerTuple
from scipy.special import jv as bessel

import data_manager

initial_state = "BCDW"
psi_initial="local"
hopping = 1.
num_sites = 2
num_electrons = num_sites//2
save_image = True

# number of ehrenfest trajectories (randomly chosen)
# use num_samples = 0 for single trajectory starting in rest on all sites
num_samples = 1

t_max = 20

omega = 0.1
gamma = [0.1, 0.4]

sim = [data_manager.Simulation(omega, gamma[i], t_max,
                               num_samples=num_samples, num_sites=num_sites,
                               symmetric_y_init=False,
                               num_electrons=num_electrons,
                               printProgress=True,
                               initial_state=initial_state,
                               psi_initial=psi_initial,
                               calcYSquared=True,
                               fullStorage=False,
                               hopping=hopping) for i in range(len(gamma))]

ref = [data_manager.RefData(initial_state, psi_initial, num_sites, omega, gamma[i]) for i in range(len(gamma))]

from matplotlib import gridspec

fig = plt.figure(figsize=(5.5, 3.))
gs = gridspec.GridSpec(2, len(gamma), fig, hspace=0, wspace=.23,
                       left=0.115, bottom=0.136, right=.985, top=.982)

#fig.suptitle("         L=%i, "%num_sites + initial_state + ", $\\omega=%gt_0$"%omega)

linestyle_ref = "dotted"

text_label = ["(a)", "(c)", "(b)", "(d)"]

# the first subplot
for i in range(len(gamma)):
	ax0 = plt.subplot(gs[0, i])
	ax1 = plt.subplot(gs[1, i], sharex=ax0)
	if i == len(gamma)-1:
		ax = ax0 if initial_state == "DCDW" else ax1
		ax.plot([], [], "k")
		ax.plot([], [], "k", linestyle=linestyle_ref)
		ax.legend(("Ehrenfest", "ED"))
	
	ax0.plot(sim[i].time(), sim[i].O_CDW())
	ax0.plot(ref[i].time(), ref[i].O_CDW(), linestyle=linestyle_ref)
	
	#ax0.set_title("$\\gamma=%gt_0$"%(sim[i].gamma))
	ax0.text(0.9, 0.12, text_label[2*i], transform=ax0.transAxes)
	ax1.text(0.9, 0.12, text_label[2*i+1], transform=ax1.transAxes)
	
	nph_ehren = ax1.plot(sim[i].time(), sim[i].nph()/(sim[i].gamma/omega)**2)
	nph_ref = ax1.plot(ref[i].time(), ref[i].nph()/(sim[i].gamma/omega)**2, linestyle="dotted")

	plt.setp(ax0.get_xticklabels(), visible=False)
	# remove last tick label for the second subplot
	#yticks = ax0.yaxis.get_major_ticks()
	#yticks[0].label1.set_visible(False)
	
	ax0.set_xlim((0, t_max))
	ax0.set_ylim((-1.02, 1.02))
	ax1.set_ylim((-0.02, 1.1))
	
	ax1.set_xlabel("$t\\cdot t_0$")
	
	if i == 0:
		ax0.set_ylabel("$\\mathcal{O}_{\\mathsf{CDW}}$")
		ax1.set_ylabel("$N_{ph}\\cdot \\omega^2/\\gamma^2$")
		#ax1.legend([(nph_ehren[0], nph_ref[0]), (nph_ehren[1], nph_ref[1])], ["occupied", "empty"], numpoints=1,
		#           handler_map={tuple: HandlerTuple(ndivide=None)})
	

if save_image:
	fig.savefig("../report_research_and_dev/figures/L=2_%s_multi_gamma_OCDW_Nph.pdf"%initial_state)
else:
	plt.show()
