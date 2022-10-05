import numpy as np
import tt

def initial_tt_purestate(n, ranks, index_initial_state):
    """
    n=[n_0,n_1,...,n_{d-1}]
    ranks=[r_0,r_1,r_2,...,r_{d-1},r_d]
    return initial state of the system density matrices in TT format
    """

    rv=[]
    for k in range(len(n)):
        mat_temp = np.zeros((ranks[k],n[k],ranks[k+1]), dtype = complex)
        mat_temp[0, index_initial_state[k], 0] = 1. 
        rv.append(mat_temp)

    rv=tt.vector.from_list(rv)
    return rv


def initial_tt_rho(n, ranks, index_initial_state_env, n_elec, vib_basis):
    """
    n=[n_0,n_1,...,n_{d-1}]
    ranks=[r_0,r_1,r_2,...,r_{d-1},r_d]
    return initial state of the system density matrices in TT format
    """

    rv=[]
    length = len(n)
    length_env = len(index_initial_state_env)
    length_elec = length_env+n_elec
    for k in range(length_env):
        mat_temp = np.zeros((ranks[k],n[k],ranks[k+1]), dtype = complex)
        mat_temp[0, index_initial_state_env[k], 0] = 1. 
        rv.append(mat_temp)
    for k in range(length_env, length_elec):
        mat_temp = np.zeros((ranks[k],n[k],ranks[k+1]), dtype = complex)
        mat_temp[0,0,0]=0.5
        mat_temp[0,3,0]=0.5
        rv.append(mat_temp)
    for k in range(length_elec, length):
        mat_temp = np.zeros((ranks[k],n[k],ranks[k+1]), dtype = complex)
        i_vib = k - length_elec
        for i_state in range(vib_basis[i_vib]):
            i_pop = (vib_basis[i_vib]+1)*i_state
            mat_temp[0,i_pop,0]=1.0/vib_basis[i_vib]
        rv.append(mat_temp)

    rv=tt.vector.from_list(rv)
    return rv


    


    


