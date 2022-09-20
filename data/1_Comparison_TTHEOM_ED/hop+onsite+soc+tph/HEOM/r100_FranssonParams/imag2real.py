import numpy as np
import itertools
import tt

def modes_imag(number_el_state, number_spin, number_lead, number_pole, head_tail_coup):
        """
        Generating the index_mode: mode in {[0:number_el_state,0:number_lead,0:number_pole+1,0:sign_pm, sin/cos]}
        i.e., sin/cos=(0,1)[s*cos+(1-s)*sin], number_lead=2(left:0, right:1),    number_el_state=1,  number_pole=1,  sign_pm=2(0:+, 1:-)
        """
        el_states=np.arange(number_el_state)
        spins=np.arange(number_spin)
        leads=np.arange(number_lead)
        poles=np.arange(number_pole)
        sign=np.arange(2)
        sin_cos=np.arange(2) #[0,1] 0:sin, 1:cos
    
        modes_imag=[]
        if head_tail_coup:
            for elem in itertools.product(spins, leads, poles, sign, sin_cos):
                modes_imag.append(list(elem))
        else:
            for elem in itertools.product(el_states, spins, leads, poles, sign, sin_cos):
                modes_imag.append(list(elem))
    
        return np.array(modes_imag)

def transform(tensor, N_env, modes, coefficients, partition_perDOF, head_tail_coup):
    """Get element of the TT-vector.
    :param index: array_like (it supports slicing).
    :returns: number -- an element of the tensor or a new tensor.
    Examples:
    Suppose that a is a 3-dimensional tt.vector of size 4 x 5 x 6
    a[1, 2, 3] returns the element with index (1, 2, 3)
    a[1, :, 1:3] returns a 2-dimensional tt.vector of size 5 x 2
    """
    d = tensor.d
    answ_cores = []
    for i in range(0, N_env, 2):
        # r0, n, r1 = cores[i].shape
        cur_core_1 = tensor.core[tensor.ps[i] - 1 : tensor.ps[i + 1] - 1]
        cur_core_1 = cur_core_1.reshape( (tensor.r[i], tensor.n[i], tensor.r[i + 1]), order='F')
        cur_core_2 = tensor.core[tensor.ps[i + 1] - 1 : tensor.ps[i + 2] - 1]
        cur_core_2 = cur_core_2.reshape( (tensor.r[i+1], tensor.n[i+1], tensor.r[i + 2]), order='F')
        cur_core = np.zeros( (tensor.r[i], 2, tensor.r[i+2]), dtype=complex)

        if head_tail_coup:
            i_spin = modes[i,0]
            i_lead = modes[i,1]
            i_pole = modes[i,2]
        else:
            i_elec = modes[i,0]
            i_spin = modes[i,1]
            i_lead = modes[i,2]
            i_pole = modes[i,3]
        
        cur_core_a = cur_core_1[:, 0, :].reshape((tensor.r[i], tensor.r[i+1]), order='F')
        cur_core_b = cur_core_2[:, 0, :].reshape((tensor.r[i+1], tensor.r[i+2]), order='F')
        cur_core[:, 0, :] = np.dot(cur_core_a, cur_core_b)        
        cur_core_a = cur_core_1[:, 1, :].reshape((tensor.r[i], tensor.r[i+1]), order='F')
        cur_core_b = cur_core_2[:, 0, :].reshape((tensor.r[i+1], tensor.r[i+2]), order='F')
        cur_core[:, 1, :] = -coefficients[i_lead][i_spin][i_pole].conjugate()*np.dot(cur_core_a, cur_core_b)       
        cur_core_a = cur_core_1[:, 0, :].reshape((tensor.r[i], tensor.r[i+1]), order='F')
        cur_core_b = cur_core_2[:, 1, :].reshape((tensor.r[i+1], tensor.r[i+2]), order='F')
        cur_core[:, 1, :] += -1j*coefficients[i_lead][i_spin][i_pole].conjugate()*np.dot(cur_core_a, cur_core_b)       
        answ_cores.append(cur_core)


    for i in range(N_env, d):
        cur_core = tensor.core[tensor.ps[i] - 1 : tensor.ps[i + 1] - 1]
        cur_core = cur_core.reshape( (tensor.r[i], tensor.n[i], tensor.r[i + 1]), order='F')
        answ_cores.append(cur_core/partition_perDOF)

    return tensor.from_list(answ_cores)

def transform_reduced(tensor, N_env, partition_perDOF):
    """Get element of the TT-vector.
    :param index: array_like (it supports slicing).
    :returns: number -- an element of the tensor or a new tensor.
    Examples:
    Suppose that a is a 3-dimensional tt.vector of size 4 x 5 x 6
    a[1, 2, 3] returns the element with index (1, 2, 3)
    a[1, :, 1:3] returns a 2-dimensional tt.vector of size 5 x 2
    """
    d = tensor.d
    answ_cores = []
    for i in range(0, N_env, 2):
        # r0, n, r1 = cores[i].shape
        cur_core_1 = tensor.core[tensor.ps[i] - 1 : tensor.ps[i + 1] - 1]
        cur_core_1 = cur_core_1.reshape( (tensor.r[i], tensor.n[i], tensor.r[i + 1]), order='F')
        cur_core_2 = tensor.core[tensor.ps[i + 1] - 1 : tensor.ps[i + 2] - 1]
        cur_core_2 = cur_core_2.reshape( (tensor.r[i+1], tensor.n[i+1], tensor.r[i + 2]), order='F')
        cur_core = np.zeros( (tensor.r[i], 2, tensor.r[i+2]), dtype=complex)

        cur_core_a = cur_core_1[:, 0, :].reshape((tensor.r[i], tensor.r[i+1]), order='F')
        cur_core_b = cur_core_2[:, 0, :].reshape((tensor.r[i+1], tensor.r[i+2]), order='F')
        cur_core[:, 0, :] = np.dot(cur_core_a, cur_core_b)        
        answ_cores.append(cur_core)


    for i in range(N_env, d):
        cur_core = tensor.core[tensor.ps[i] - 1 : tensor.ps[i + 1] - 1]
        cur_core = cur_core.reshape( (tensor.r[i], tensor.n[i], tensor.r[i + 1]), order='F')
        answ_cores.append(cur_core/partition_perDOF)

    return tensor.from_list(answ_cores)

