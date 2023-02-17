import numpy as np
import tt

def reduced_tensor(tensor,index_1, index_2):
    """Get element of the TT-vector.
    :param index: array_like (it supports slicing).
    :returns: number -- an element of the tensor or a new tensor.
    Examples:
    Suppose that a is a 3-dimensional tt.vector of size 4 x 5 x 6
    a[1, 2, 3] returns the element with index (1, 2, 3)
    a[1, :, 1:3] returns a 2-dimensional tt.vector of size 5 x 2
    """
    d = tensor.d
    if (len(index_1) + len(index_2))  > d:
        print("Incorrect index length.")
        return
    running_fact = None
    for i in range(len(index_1)):
        # r0, n, r1 = cores[i].shape
        cur_core = tensor.core[tensor.ps[i] - 1 : tensor.ps[i + 1] - 1]
        cur_core = cur_core.reshape( (tensor.r[i], tensor.n[i], tensor.r[i + 1]), order='F')
        cur_core = cur_core[:, index_1[i], :].reshape((tensor.r[i], -1), order='F')
        if running_fact is not None:
            cur_core = np.dot(running_fact, cur_core)
        running_fact  = cur_core.reshape((1, tensor.r[i + 1]), order='F')

    answ_cores = []
    for i in range(len(index_1), d - len(index_2)):
        cur_core = tensor.core[tensor.ps[i] - 1 : tensor.ps[i + 1] - 1]
        cur_core = cur_core.reshape( (tensor.r[i], tensor.n[i], tensor.r[i + 1]), order='F')
        answ_cores.append(cur_core)

    answ_cores[0] = np.tensordot(running_fact, answ_cores[0], (1,0))

    running_fact = None
    for i in range(d-1, d-1 - len(index_2), -1):
        # r0, n, r1 = cores[i].shape
        cur_core = tensor.core[tensor.ps[i] - 1 : tensor.ps[i + 1] - 1]
        cur_core = cur_core.reshape( (tensor.r[i], tensor.n[i], tensor.r[i + 1]), order='F')
        i_site=len(index_2)-d+i
        cur_core = cur_core[:, index_2[i_site], :].reshape((-1, tensor.r[i+1]), order='F')
        if running_fact is not None:
            cur_core = np.dot(cur_core, running_fact)
        running_fact  = cur_core.reshape((tensor.r[i], 1), order='F')

    answ_cores[-1] = np.tensordot(answ_cores[-1], running_fact, (2,0))

    return tensor.from_list(answ_cores)

def reduced_tensor_L(tensor,index):
    """Get element of the TT-vector.
    :param index: array_like (it supports slicing).
    :returns: number -- an element of the tensor or a new tensor.
    Examples:
    Suppose that a is a 3-dimensional tt.vector of size 4 x 5 x 6
    a[1, 2, 3] returns the element with index (1, 2, 3)
    a[1, :, 1:3] returns a 2-dimensional tt.vector of size 5 x 2
    """
    d = tensor.d
    if len(index)  > d:
        print("Incorrect index length.")
        return
    running_fact = None
    answ_cores = []
    for i in range(len(index)):
        # r0, n, r1 = cores[i].shape
        cur_core = tensor.core[tensor.ps[i] - 1 : tensor.ps[i + 1] - 1]
        cur_core = cur_core.reshape( (tensor.r[i], tensor.n[i], tensor.r[i + 1]), order='F')
        cur_core = cur_core[:, index[i], :].reshape((tensor.r[i], -1), order='F')
        if running_fact is not None:
            cur_core = np.dot(running_fact, cur_core)
        running_fact  = cur_core.reshape((1, tensor.r[i + 1]), order='F')

    for i in range(len(index), d):
        cur_core = tensor.core[tensor.ps[i] - 1 : tensor.ps[i + 1] - 1]
        cur_core = cur_core.reshape( (tensor.r[i], tensor.n[i], tensor.r[i + 1]), order='F')
        answ_cores.append(cur_core)

    answ_cores[0] = np.tensordot(running_fact, answ_cores[0], (1,0))

    return tensor.from_list(answ_cores)


def reduced_tensor_R(tensor,index):
    """Get element of the TT-vector.
    :param index: array_like (it supports slicing).
    :returns: number -- an element of the tensor or a new tensor.
    Examples:
    Suppose that a is a 3-dimensional tt.vector of size 4 x 5 x 6
    a[1, 2, 3] returns the element with index (1, 2, 3)
    a[1, :, 1:3] returns a 2-dimensional tt.vector of size 5 x 2
    """
    d = tensor.d
    if len(index)  > d:
        print("Incorrect index length.")
        return

    answ_cores = []
    for i in range(d - len(index)):
        cur_core = tensor.core[tensor.ps[i] - 1 : tensor.ps[i + 1] - 1]
        cur_core = cur_core.reshape( (tensor.r[i], tensor.n[i], tensor.r[i + 1]), order='F')
        answ_cores.append(cur_core)

    running_fact = None
    for i in range(d-1, d-1 - len(index), -1):
        # r0, n, r1 = cores[i].shape
        cur_core = tensor.core[tensor.ps[i] - 1 : tensor.ps[i + 1] - 1]
        cur_core = cur_core.reshape( (tensor.r[i], tensor.n[i], tensor.r[i + 1]), order='F')
        i_site=len(index)-d+i
        cur_core = cur_core[:, index[i_site], :].reshape((-1, tensor.r[i+1]), order='F')
        if running_fact is not None:
            cur_core = np.dot(cur_core, running_fact)
        running_fact  = cur_core.reshape((tensor.r[i], 1), order='F')

    answ_cores[-1] = np.tensordot(answ_cores[-1], running_fact, (2,0))

    return tensor.from_list(answ_cores)


