import numpy as np
import tt

def OP1D_tt(mat, I_mat, i, dim):
    """
    return I(1) \otimes I(2) \otimes ...  \otimes I(i-1) \otimes mat  \otimes  I(i+1) \otimes  I(dim)
    """
    W=mat
    for j in range(i,0,-1):
        W = tt.kron( I_mat[j-1], W )
    for j in range(i+1, dim):
        W = tt.kron( W, I_mat[j] )
    return W

def OP1D_tt_fermi(mat, I_mat_L, I_mat_R, i, dim):
    """
    return I_L(1) \otimes I_L(2) \otimes ... \otimes I_L(i-1) \otimes mat \otimes I_R(i+1) \otimes I_R(dim)
    """
    W=mat
    for j in range(i,0,-1):
        W = tt.kron( I_mat_L[j-1], W )
    for j in range(i+1, dim):
        W = tt.kron( W, I_mat_R[j] )
    return W




