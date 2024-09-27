from constants import *
from input_para_vib import *
import numpy as np
import Pade, index
import tt, tt.ksl
import mpo_heom, init_mps
import op, gettensor, imag2real
from time import perf_counter, asctime
import os, sys, re
import scipy.io
from tempfile import TemporaryFile


rank_max = 120 

filename_appendix='_Hop+Onsite_bias_{V:.2f}_dt_{dt}_rank_{rank_max}'.format(V=Phi,dt=dt*au2fs,rank_max=rank_max)
#filename_appendix='_Hop+Onsite_bias_'+str(Phi*au2eV)+'_dt_'+str(dt*au2fs)+'_rank_'+str(rank_max)
output_info = open( 'simulation_info'+filename_appendix+'.txt', 'w')
t_start_program = asctime();    t_end = perf_counter(); output_info.write( 'The simulation starts at '+str(t_start_program)+'\n')

pade_fermi = Pade.Pade_decomposition(N_fermipoles, 'Fermi')
coefficients_fermi, exponents_fermi=pade_fermi.get_coefficients_exponents(temperature_fermibath, broadenings_fermibath, N_leads, N_spin, Chem_potentials,  WBL)
exponents_imag=exponents_fermi[:,:,:,0].real
coefficients_imag=(1j*coefficients_fermi).real
for i_lead in range(N_leads):
    for i_spin in range(N_spin):
        coefficients_imag[i_lead, i_spin, 0]=broadenings_fermibath[i_lead][i_spin]/4./np.cos(broadenings_fermibath[i_lead][i_spin]/(2.*temperature_fermibath[i_lead][i_spin]))

if not WBL: N_fermipoles += 1
if head_tail_coup:
    N_env = N_spin * N_leads * N_fermipoles * 2 + n_dummy
else:
    N_env = Dim_elecDOF * N_spin * N_leads * N_fermipoles * 2
index_env0 = np.zeros(N_env, dtype=int)
index_env0_imagtime = np.zeros(N_env*2, dtype=int)

hierarchy=index.Hierarchy_index(Dim_elecDOF, N_spin, N_leads, N_fermipoles, head_tail_coup, 2)
indices_realtime=hierarchy.indices_real_generation()
indices_imagtime=hierarchy.indices_imag_generation()

pop_el_output = re.search(r"\[([a-z0-9_]+)\]", Outputs.lower())
if ( pop_el_output != None ): pop_el_output = pop_el_output.group(0)[1:-1] 
el_states_output = []
if pop_el_output == 'all':
    dim_fockspace = 2**(Dim_elecDOF*N_spin)
elif ( pop_el_output != None ):
    el_states_output = np.loadtxt(pop_el_output+'.out', dtype=int)
    el_states_output = el_states_output.reshape((-1, Dim_elecDOF*N_spin))
    dim_fockspace=len(el_states_output)

pop_vib_output = re.search(r"\{([a-z0-9_]+)\}", Outputs.lower())
if ( pop_vib_output != None ): pop_vib_output = pop_vib_output.group(0)[1:-1] 
vib_states_output = []
if pop_vib_output == 'all':
    dim_vibfockspace = np.prod(vibbasis)
elif ( pop_vib_output != None ):
    vib_states_output = np.loadtxt(pop_vib_output+'.out', dtype=int)
    vib_states_output = vib_states_output.reshape((-1, Dim_vibDOF))
    dim_vibfockspace=len(vib_states_output)

n_vib_output = False
if 'n_vib' in Outputs.lower(): n_vib_output = True

current_output = False
if 'current' in Outputs.lower(): current_output = True

if Init_propagation:
    if ('ground_state' in Init_state.lower()) or ('initial_state' in Init_state.lower()):
        mpo = mpo_heom.TTOP(Dim_elecDOF, N_spin, Dim_vibDOF, vibbasis, N_leads, N_fermipoles, n_dummy, beta, SPIN_ALT_ORDER=Spin_alt_order, HEAD_TAIL_COUP = head_tail_coup, TWIN_SPACE = True, eps=eps)
    elif 'equilibrium' in Init_state.lower():
        mpo = mpo_heom.TTOP(Dim_elecDOF, N_spin, Dim_vibDOF, vibbasis, N_leads, N_fermipoles, n_dummy, beta, SPIN_ALT_ORDER=Spin_alt_order, HEAD_TAIL_COUP = head_tail_coup, TWIN_SPACE = False, eps=eps)

    #mpo.OP_hamiltonian_SIAM(onsite_energy, Couloumb_energy, exchange_energy, vib_freq, vib_elec_coup) 
    mpo.OP_hamiltonian_CISS(onsite_energy, Couloumb_energy, exchange_energy, spin_orbit_coupling, nnhopping_coupling_ph, spin_orbit_coupling_ph, vib_freq, vib_elec_coup)
    II, Pop_el, Pop_vib, N_vib, N_site = mpo.OP_observables(EL_STATES_OUTPUT = el_states_output, VIB_STATES_OUTPUT = vib_states_output, POP_EL_OUTPUT= pop_el_output, POP_VIB_OUTPUT= pop_vib_output, N_VIB_OUTPUT= n_vib_output)

    XOP_realtime, Current = mpo.HQME_fermi_env_realtime(coefficients_fermi, exponents_fermi, elec_lead_couplings, WBL=WBL)
    n_realtime=XOP_realtime.n  # n_realtime is the overall size of the basis set. n_realtime[0] contains the depth of the hierarchy.
    dof_realtime = XOP_realtime.tt.d # dof_realtime contains the number of physical variables and the number of variables included in the tilde space.

    XOP_list = tt.matrix.to_list(XOP_realtime)
    #np.savez("npz/XOP"+filename_appendix+".npz",*XOP_list)
    XOP_list_tt=tt.vector.to_list(XOP_realtime.tt)
    tpz = TemporaryFile()
    np.savez(tpz,*XOP_list_tt)
    _ = tpz.seek(0)
    scipy.io.savemat("../tamen/XOP.mat",np.load(tpz))
    
    II_list = tt.vector.to_list(II)
    #np.savez("npz/II"+filename_appendix+".npz",*II_list)

    if ('ground_state' in Init_state.lower()) or ('initial_state' in Init_state.lower()):

        ranks_realtime = np.zeros(dof_realtime+1, dtype=int)
        ranks_realtime[0] = 1
        ranks_realtime[dof_realtime] = 1
        for i in range(1, dof_realtime):
            ranks_realtime[i] = min(rank_max, ranks_realtime[i-1]*n_realtime[i-1])
        for i in range(dof_realtime-1, 0, -1):
            ranks_realtime[i] = min(ranks_realtime[i], ranks_realtime[i+1]*n_realtime[i])
        output_info.write('ranks_realtime:\t'+str(ranks_realtime)+'\n')
        if 'ground_state' in Init_state.lower():
            initial_state = np.zeros(len(n_realtime), dtype=int)
        else:
            initial_state_system = np.loadtxt(Init_state+".in")
            initial_state = np.zeros(len(n_realtime), dtype=int)
            initial_state[N_env:N_env+2*(Dim_elecDOF*N_spin+Dim_vibDOF)] = np.array( [ [i,j] for i,j in zip(initial_state_system, initial_state_system) ] ).ravel()
            
        rho0=init_mps.initial_tt_purestate(n_realtime, ranks_realtime, initial_state)
    elif 'equilibrium' in Init_state.lower():

        XOP_imagtime = mpo.HQME_fermi_env_imagtime(coefficients_imag, exponents_imag, elec_lead_couplings, beta)
        n_imagtime = XOP_imagtime.n
        dof_imagtime = XOP_imagtime.tt.d

        ranks_imagtime = np.zeros(dof_imagtime+1, dtype=int)
        ranks_imagtime[0] = 1
        ranks_imagtime[dof_imagtime] = 1
        for i in range(1, dof_imagtime):
            ranks_imagtime[i] = min(rank_max, ranks_imagtime[i-1]*n_imagtime[i-1])
        for i in range(dof_realtime-1, 0, -1):
            ranks_imagtime[i] = min(ranks_imagtime[i], ranks_imagtime[i+1]*n_imagtime[i])
        output_info.write('ranks_imagtime:\t'+str(ranks_imagtime)+'\n')

        initial_state_env = np.zeros((N_env*2), dtype=int)
        rho_imag = init_mps.initial_tt_rho(n_imagtime, ranks_imagtime, initial_state_env, Dim_elecDOF*N_spin, vibbasis)

        xtau = np.arange(1, 1+Ntaus, dtype=float) * dtau 
        for i_tau in range(Ntaus):
            XOP_imagtime = mpo.HQME_fermi_env_imagtime(coefficients_imag, exponents_imag, elec_lead_couplings, xtau[i_tau])
            rho_imag = tt.ksl.ksl(XOP_imagtime, rho_imag, dtau)

        #f_adm_r = open('adm_imag_r'+filename_appendix+'.txt', 'w')
        #f_adm_i = open('adm_imag_i'+filename_appendix+'.txt', 'w')
        #for i in range(len(indices_imagtime)):
        #    index_env = list(map(int,indices_imagtime[i]))
        #    rho_sys = gettensor.reduced_tensor_L(rho_imag, index_env)
        #    f_adm_r.write( str(i) + '\t' )
        #    f_adm_i.write( str(i) + '\t' )
        #    for i_elecstate in range(dim_fockspace**2):
        #        pop = tt.dot(Pop_el[i_elecstate], rho_sys)
        #        f_adm_r.write( str(pop.real) + '\t' )
        #        f_adm_i.write( str(pop.imag) + '\t' )
        #    f_adm_r.write( '\n' )
        #    f_adm_i.write( '\n' )

        rho_sys = gettensor.reduced_tensor_L(rho_imag, index_env0_imagtime)
        norm = tt.dot(mpo.I1_sys, rho_sys) 
        partition = norm.real
        if 'reduced' in Init_state.lower():
            rho0 = imag2real.transform_reduced(rho_imag, N_env*2, partition**(1./(Dim_elecDOF*N_spin+Dim_vibDOF)))
        else:
            modes = imag2real.modes_imag(Dim_elecDOF, N_spin, N_leads, N_fermipoles, head_tail_coup)
            rho0 = imag2real.transform(rho_imag, N_env*2, modes, coefficients_fermi, partition**(1./(Dim_elecDOF*N_spin+Dim_vibDOF)), head_tail_coup)

        #f_adm_r = open('adm_tran_r'+filename_appendix+'.txt', 'w')
        #f_adm_i = open('adm_tran_i'+filename_appendix+'.txt', 'w')
        #for i in range(len(indices_realtime)):
        #    index_env = list(map(int,indices_realtime[i]))
        #    rho_sys = gettensor.reduced_tensor_L(rho0, index_env)
        #    f_adm_r.write( str(i) + '\t' )
        #    f_adm_i.write( str(i) + '\t' )
        #    for i_elecstate in range(dim_fockspace**2):
        #        pop = tt.dot(Pop_el[i_elecstate], rho_sys)
        #        f_adm_r.write( str(pop.real) + '\t' )
        #        f_adm_i.write( str(pop.imag) + '\t' )
        #    f_adm_r.write( '\n' )
        #    f_adm_i.write( '\n' )

    rho_sys = gettensor.reduced_tensor_L(rho0, index_env0)
    f_norm = open('norm_realtime'+filename_appendix+'.txt', 'w')
    f_norm.write( str(0.) + '\t' )
    norm = tt.dot(II, rho_sys) 
    f_norm.write(str(norm.real) + '\t' + str(norm.imag) + '\n')
    f_norm.close()

    if pop_el_output != None:
        f_pop_elec = open('pop_elec_realtime'+filename_appendix+'.txt', 'w')
        f_pop_elec.write( str(0.) + '\t' )
        for i_elecstate in range(dim_fockspace):
            pop = tt.dot(Pop_el[i_elecstate], rho_sys)
            f_pop_elec.write( str(pop.real) + '\t' )
            Pop_el_list = tt.vector.to_list(Pop_el[i_elecstate])
            #np.savez("npz/Pop_el_"+str(i_elecstate)+filename_appendix+".npz",*Pop_el_list)
            #Pop = np.load("npz/Pop_el_"+str(i_elecstate)+filename_appendix+".npz")
            #scipy.io.savemat("npz/Pop_el_"+str(i_elecstate)+filename_appendix+".mat",Pop)

        f_pop_elec.write('\n')
        f_pop_elec.close()
    if pop_vib_output != None:
        f_pop_vib = open('pop_vib_realtime'+filename_appendix+'.txt', 'w')
        f_pop_vib.write( str(0.) + '\t' )
        for i_vibstate in range(dim_vibfockspace):
            pop = tt.dot(Pop_vib[i_vibstate], rho_sys)
            f_pop_vib.write( str(pop.real) + '\t' )
            Pop_vib_list = tt.vector.to_list(Pop_vib[i_vibstate])
            #np.savez("npz/Pop_vib_"+str(i_vibstate)+filename_appendix+".npz",*Pop_vib_list)
            #Pop = np.load("npz/Pop_vib_"+str(i_vibstate)+filename_appendix+".npz")
            #scipy.io.savemat("npz/Pop_vib_"+str(i_vibstate)+filename_appendix+".mat",Pop)
        f_pop_vib.write('\n')
        f_pop_vib.close()
    if n_vib_output:
        f_n_vib = open('n_vib_realtime'+filename_appendix+'.txt', 'w')
        f_n_vib.write( str(0.) + '\t' )
        for i_vib in range(Dim_vibDOF):
            nvib= tt.dot(N_vib[i_vib], rho_sys)
            f_n_vib.write( str(nvib.real) + '\t' )
            N_vib_list = tt.vector.to_list(N_vib[i_vib])
            #np.savez("npz/N_vib_"+str(i_vib)+filename_appendix+".npz",*N_vib_list)
            #Nvib = np.load("npz/N_vib_"+str(i_vib)+filename_appendix+".npz")
            #scipy.io.savemat("npz/N_vib_"+str(i_vib)+filename_appendix+".mat",Nvib)
        f_n_vib.write('\n')
        f_n_vib.close()
    f_n_site = open('n_site_realtime'+filename_appendix+'.txt', 'w')
    f_n_site.write( str(0.) + '\t' )
    for i_elec in range(Dim_elecDOF*N_spin):
        nsite= tt.dot(N_site[i_elec], rho_sys)
        f_n_site.write( str(nsite.real) + '\t' )
        N_site_list = tt.vector.to_list(N_site[i_elec])
        #np.savez("npz/N_site_"+str(i_elec)+filename_appendix+".npz",*N_site_list)
        #Nsite = np.load("npz/N_site_"+str(i_elec)+filename_appendix+".npz")
        #scipy.io.savemat("npz/N_site_"+str(i_elec)+filename_appendix+".mat",Nsite)

    f_n_site.write('\n')
    f_n_site.close()
    if current_output:
        f_current = open('current_realtime'+filename_appendix+'.txt', 'w')
        f_current.write("#(0)time, currents:\t")
        for i_spin in range(N_spin):
            for i_lead in range(N_leads):
                i_state = i_spin*N_leads+i_lead
                f_current.write("("+str(i_state+1)+") lead_"+str(i_lead)+"_spin_"+str(i_spin)+"\t")
        f_current.write("\n")
        f_current.write( str(0.) + '\t' )
        current=np.zeros((N_leads*N_spin), dtype=float)
        for i_spin in range(N_spin):
            for i_lead in range(N_leads):
                i_state = i_spin*N_leads+i_lead
                Current_list = tt.vector.to_list(Current[i_state])
                np.savez("npz/Current_"+str(i_spin)+str(i_lead)+filename_appendix+".npz",*Current_list)
                current[i_state] = np.real( tt.dot(Current[i_state], rho0))
                f_current.write( str( current[i_state] * au2ampere ) +'\t' )
                Cur = np.load("npz/Current_"+str(i_spin)+str(i_lead)+filename_appendix+".npz")
                scipy.io.savemat("../tamen/CurrOp"+str(i_spin)+str(i_lead)+".mat",Cur)
        f_current.write('\n')
        f_current.close()

else:
    XOP_list = np.load("npz/XOP"+filename_appendix+".npz", mmap_mode='r')
    data =[XOP_list[i] for i in XOP_list]
    XOP_realtime = tt.matrix.from_list(data)
    II_list = np.load("npz/II"+filename_appendix+".npz")
    data =[II_list[i] for i in II_list]
    II = tt.vector.from_list(data)

    n_realtime=XOP_realtime.n  # n_realtime is the overall size of the basis set. n_realtime[0] contains the depth of the hierarchy.
    dof_realtime = XOP_realtime.tt.d # dof_realtime contains the number of physical variables and the number of variables included in the tilde space.

    rho0_list = np.load("npz/rho0"+filename_appendix+".npz")
    data =[rho0_list[i] for i in rho0_list]
    rho0 = tt.vector.from_list(data)
    ranks_realtime = rho0.r

    if pop_el_output != None:
        Pop_el = [None] * dim_fockspace
        for i_elecstate in range(dim_fockspace):
            Pop_el_list = np.load("npz/Pop_el_"+str(i_elecstate)+filename_appendix+".npz")
            
            data =[Pop_el_list[i] for i in Pop_el_list]
            Pop_el[i_elecstate] = tt.vector.from_list(data)
    if pop_vib_output != None:
        Pop_vib = [None] * dim_vibfockspace
        for i_vibstate in range(dim_vibfockspace):
            Pop_vib_list = np.load("npz/Pop_vib_"+str(i_vibstate)+filename_appendix+".npz")
            data =[Pop_vib_list[i] for i in Pop_vib_list]
            Pop_vib[i_vibstate] = tt.vector.from_list(data)
    if n_vib_output:
        N_vib = [None] * Dim_vibDOF
        for i_vib in range(Dim_vibDOF):
            N_vib_list = np.load("npz/N_vib_"+str(i_vib)+filename_appendix+".npz")
            data =[N_vib_list[i] for i in N_vib_list]
            N_vib[i_vib] = tt.vector.from_list(data)
    N_site = [None] * Dim_elecDOF * N_spin
    for i_elec in range(Dim_elecDOF * N_spin):
        N_site_list = np.load("npz/N_site_"+str(i_elec)+filename_appendix+".npz")
        data =[N_site_list[i] for i in N_site_list]
        N_site[i_elec] = tt.vector.from_list(data)
    if current_output:
        Current = [None] * N_spin * N_leads
        for i_spin in range(N_spin):
            for i_lead in range(N_leads):
                i_state = i_spin*N_leads+i_lead
                Current_list = np.load("npz/Current_"+str(i_spin)+str(i_lead)+filename_appendix+".npz")
                data = [Current_list[i] for i in Current_list]
                Current[i_state] = tt.vector.from_list(data)

output_info.write( 'n_realtime is the overall size of the basis set:\t'+ str(n_realtime) +'\n')
output_info.write( 'dof_realtime contains the number of physical variables:\t'+ str(dof_realtime)+'\n')
output_info.write( 'XOP_realtime: '+str(XOP_realtime) +'\n')
output_info.write( 'XOP_realtime.norm:\t'+ str(XOP_realtime.norm()) +'\n')
output_info.write( 'II: '+str(II) +'\n')
output_info.write( 'II.norm:\t'+ str(II.norm()) +'\n')
t_start = t_end;    t_end = perf_counter(); output_info.write("Elapsed time of tensor generation: " + str(t_end-t_start) +'\n')

xtime = np.arange(1,1+Nsteps,dtype=float) * dt * au2fs + t_in * au2fs

time_section = [Nsteps//Nsections * i for i in range(Nsections)] 
time_section.append(Nsteps)

rho0_list = tt.vector.to_list(rho0)
np.savez("npz/rho0"+filename_appendix+".npz",*rho0_list)
r0 = np.load("npz/rho0"+filename_appendix+".npz")
scipy.io.savemat("../tamen/RHO.mat",r0)

for i_section in range(Nsections):
    f_norm = open('norm_realtime'+filename_appendix+'.txt', 'a+')
    f_n_site = open('n_site_realtime'+filename_appendix+'.txt', 'a+')
    if pop_el_output != None: f_pop_elec = open('pop_elec_realtime'+filename_appendix+'.txt', 'a+')
    if pop_vib_output != None: f_pop_vib = open('pop_vib_realtime'+filename_appendix+'.txt', 'a+')
    if n_vib_output: f_n_vib = open('n_vib_realtime'+filename_appendix+'.txt', 'a+')
    if current_output: f_current = open('current_realtime'+filename_appendix+'.txt', 'a+')
    

    for i_t in range(time_section[i_section], time_section[i_section+1]): 
        rho0 = tt.ksl.ksl(XOP_realtime, rho0, dt)

        rho_sys = gettensor.reduced_tensor_L(rho0, index_env0)

        f_norm.write( str(xtime[i_t]) + '\t' )
        norm = tt.dot(II, rho_sys) 
        f_norm.write(str(norm.real) + '\t' + str(norm.imag) + '\n')
        if pop_el_output != None:
            f_pop_elec.write( str(xtime[i_t]) + '\t' )
            for i_elecstate in range(dim_fockspace):
                pop = tt.dot(Pop_el[i_elecstate], rho_sys)
                f_pop_elec.write( str(pop.real) + '\t' )
            f_pop_elec.write('\n')
        if pop_vib_output != None:
            f_pop_vib.write( str(xtime[i_t]) + '\t' )
            for i_vibstate in range(dim_vibfockspace):
                pop = tt.dot(Pop_vib[i_vibstate], rho_sys)
                f_pop_vib.write( str(pop.real) + '\t' )
            f_pop_vib.write('\n')
        if n_vib_output:
            f_n_vib.write( str(xtime[i_t]) + '\t' )
            for i_vib in range(Dim_vibDOF):
                nvib= tt.dot(N_vib[i_vib], rho_sys)
                f_n_vib.write( str(nvib.real) + '\t' )
            f_n_vib.write('\n')
        f_n_site.write( str(xtime[i_t]) + '\t' )
        for i_elec in range(Dim_elecDOF*N_spin):
            nsite= tt.dot(N_site[i_elec], rho_sys)
            f_n_site.write( str(nsite.real) + '\t' )
        f_n_site.write('\n')
        if current_output:
            f_current.write( str(xtime[i_t]) + '\t' )
            current=np.zeros((N_leads*N_spin), dtype=float)
            for i_spin in range(N_spin):
                for i_lead in range(N_leads):
                    i_state = i_spin*N_leads+i_lead
                    current[i_state] = np.real( tt.dot(Current[i_state], rho0))
                    f_current.write( str(current[i_state] * au2ampere) + '\t')
            f_current.write('\n')
    f_norm.close()
    f_n_site.close()
    if pop_el_output != None: f_pop_elec.close()
    if pop_vib_output != None: f_pop_vib.close()
    if n_vib_output: f_n_vib.close()
    if current_output: f_current.close()

#f_adm_r = open('adm_real_r'+filename_appendix+'.txt', 'w')
#f_adm_i = open('adm_real_i'+filename_appendix+'.txt', 'w')
#for i in range(len(indices_realtime)):
#    index_env = list(map(int,indices_realtime[i]))
#    rho_sys = gettensor.reduced_tensor_L(rho0, index_env)
#    f_adm_r.write( str(i) + '\t' )
#    f_adm_i.write( str(i) + '\t' )
#    for i_elecstate in range(dim_fockspace**2):
#        pop = tt.dot(Pop_el[i_elecstate], rho_sys)
#        f_adm_r.write( str(pop.real) + '\t' )
#        f_adm_i.write( str(pop.imag) + '\t' )
#    f_adm_r.write( '\n' )
#    f_adm_i.write( '\n' )
    
r0 = np.load("npz/rho0"+filename_appendix+".npz")
scipy.io.savemat("npz/rho0"+filename_appendix+".mat",r0)
t_start = t_end;    t_end = perf_counter(); output_info.write("Elapsed time of propagation: " + str(t_end-t_start) +'\n')
t_end = asctime();  output_info.write('The simulation ends at '+str(t_end)+'\n')
output_info.close()

