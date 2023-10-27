import numpy as np
import tt
import op
import itertools

class TTOP():
    def __init__(self, n_elec, n_spin, n_vib, vib_basis, n_lead, n_pole_fermi, n_dummy, beta, SPIN_ALT_ORDER = True, HEAD_TAIL_COUP = False, TWIN_SPACE = True, eps = 1.e-12):
        # eps: accuracy of tensor-train rounding
        self.n_elec = n_elec
        self.n_spin = n_spin
        self.n_vib = n_vib
        self.vib_basis = vib_basis
        self.n_lead = n_lead
        self.n_pole_fermi = n_pole_fermi
        self.n_dummy = n_dummy
        self.beta = beta
        self.SPIN_ALT_ORDER = SPIN_ALT_ORDER
        self.HEAD_TAIL_COUP = HEAD_TAIL_COUP
        self.TWIN_SPACE = TWIN_SPACE
        self.eps = eps

        self.mat_creation = tt.matrix( np.array([[0,0],[1,0]]), self.eps )
        self.mat_annihilation = tt.matrix( np.array([[0,1],[0,0]]), self.eps )
        self.mat_occupation = tt.matrix( np.array([[0,0],[0,1]]), self.eps )

        m_creation_L = np.kron( np.array([[0,0],[1,0]]), np.eye(2) )
        self.mat_creation_L = tt.matrix( m_creation_L, self.eps )
        m_creation_R = np.kron( np.eye(2), np.array([[0,0],[1,0]]) )
        self.mat_creation_R = tt.matrix( m_creation_R, self.eps )
        m_annihilation_L = np.kron( np.array([[0,1],[0,0]]), np.eye(2) )
        self.mat_annihilation_L = tt.matrix( m_annihilation_L, self.eps )
        m_annihilation_R = np.kron( np.eye(2), np.array([[0,1],[0,0]]) )
        self.mat_annihilation_R = tt.matrix( m_annihilation_R, self.eps )
        mat_occupation_L = np.kron( np.array([[0,0],[0,1]]), np.eye(2) )
        self.mat_occupation_L = tt.matrix( mat_occupation_L, self.eps )
        mat_occupation_R = np.kron( np.eye(2), np.array([[0,0],[0,1]]) )
        self.mat_occupation_R = tt.matrix( mat_occupation_R, self.eps )

        # Identity operators for system electronic DoFs in twin space
        self.I_elec_list = []
        self.I_elec_sign_L_list = []
        self.I_elec_sign_R_list = []
        if self.TWIN_SPACE:
            for i in range( 2 * self.n_elec * self.n_spin):
                self.I_elec_list.append( tt.eye( [2] ))
            for i in range(self.n_elec * self.n_spin):
                self.I_elec_sign_L_list.append( tt.matrix( np.array([[1,0],[0,-1]]), eps) )
                self.I_elec_sign_L_list.append( tt.eye( [2] ) )
                self.I_elec_sign_R_list.append( tt.eye( [2] ) )
                self.I_elec_sign_R_list.append( tt.matrix( np.array([[1,0],[0,-1]]), eps) )
            self.II_elec = tt.eye( [2] * 2 * self.n_elec * self.n_spin ) 
        else:
            for i in range(self.n_elec * self.n_spin):
                self.I_elec_list.append( tt.eye( [4] ))
            for i in range(self.n_elec * self.n_spin):
                mat = tt.matrix( np.kron(np.array([[1,0],[0,-1]]), np.eye(2)), self.eps )
                self.I_elec_sign_L_list.append( mat )
                mat = tt.matrix( np.kron(np.eye(2), np.array([[1,0],[0,-1]])), self.eps) 
                self.I_elec_sign_R_list.append( mat )
            self.II_elec = tt.eye( [4] * self.n_elec * self.n_spin ) 
    
        # Identity operators for system vibrational DoFs in twin space
        self.I_vib_list = []
        if self.TWIN_SPACE:
            for i in range(self.n_vib):
                self.I_vib_list.append( tt.eye([self.vib_basis[i]]) )
                self.I_vib_list.append( tt.eye([self.vib_basis[i]]) )
            ds_vib_basis = np.array( [ [i,j] for i,j in zip(vib_basis, vib_basis) ] ).ravel()
            self.II_vib = tt.eye( ds_vib_basis )
        else:
            for i in range(self.n_vib):
                self.I_vib_list.append( tt.eye([self.vib_basis[i]**2]) )
            ds_vib_basis = np.array( [ i**2 for i in self.vib_basis] )
            self.II_vib = tt.eye( ds_vib_basis )
        
        self.II_dummy = tt.eye([2] * 2*self.n_dummy)
        self.II_sys = tt.kron(self.II_elec, tt.kron(self.II_vib, self.II_dummy)) 

        # Unit vector for system electronic DoFs
        self.I1_elec = None
        if self.TWIN_SPACE:
            for i in range(self.n_elec * self.n_spin):
                Id = None
                for j in range(2):
                    nl =  np.zeros((2), dtype=complex) # |00> + |11>
                    nl[j] = 1.+0.0j
                    I2 = tt.vector( nl )
                    Id += tt.kron(I2, I2)
                self.I1_elec = tt.kron(self.I1_elec, Id)
        else:
            for i in range(self.n_elec * self.n_spin):
                Id = None
                nl =  np.zeros((4), dtype=complex) # |00> + |11>
                nl[0] = 1.+0.0j
                nl[3] = 1.+0.0j
                Id = tt.vector( nl )
                self.I1_elec = tt.kron(self.I1_elec, Id)
        self.I1_elec = self.I1_elec.round(self.eps)
    
        # Unit vector for system vibrational DoFs
        self.I1_vib = None
        if self.TWIN_SPACE:
            for i in range(self.n_vib):
                Id = None
                for j in range(self.vib_basis[i]):
                    nl =  np.zeros((self.vib_basis[i]), dtype=complex) # |00> + |11> + |22> + ...
                    nl[j] = 1.+0.0j
                    I2 = tt.vector( nl )  
                    Id += tt.kron(I2, I2)
                self.I1_vib = tt.kron(self.I1_vib, Id)
        else:
            for i in range(self.n_vib):
                Id = None
                nl =  np.zeros((self.vib_basis[i]**2), dtype=complex) # |00> + |11> + |22> + ...
                for j in range(self.vib_basis[i]):
                    i_state = j * self.vib_basis[i] + j
                    nl[i_state] = 1.+0.0j
                Id = tt.vector( nl )  
                self.I1_vib = tt.kron(self.I1_vib, Id)
        self.I1_vib = self.I1_vib.round(self.eps)
        self.I1_dummy = None
        for i in range(self.n_dummy):
            Id = None
            for j in range(2):
                nl =  np.zeros((2), dtype=complex) # |00> + |11>
                nl[j] = 1.+0.0j
                I2 = tt.vector( nl )
                Id += tt.kron(I2, I2)
            self.I1_dummy = tt.kron(self.I1_dummy, Id)
        
        self.I1_sys = tt.kron(self.I1_elec, tt.kron(self.I1_vib, self.I1_dummy))
        self.I1_sys = self.I1_sys.round(self.eps)
        if self.HEAD_TAIL_COUP:
            d_creation_L = []
            d_creation_R = []
            d_annihilation_L = []
            d_annihilation_R = []
            for i_spin in range(self.n_spin):
                if self.SPIN_ALT_ORDER:
                    i_state = i_spin
                else:
                    i_state = i_spin*self.n_elec
                if self.TWIN_SPACE:
                    mat = op.OP1D_tt_fermi(self.mat_creation, self.I_elec_sign_L_list, self.I_elec_list, 2*i_state, 2*self.n_elec*self.n_spin)
                else:
                    mat = op.OP1D_tt_fermi(self.mat_creation_L, self.I_elec_sign_L_list, self.I_elec_list, i_state, self.n_elec*self.n_spin)
                d_creation_L.append( tt.kron(mat, tt.kron(self.II_vib,self.II_dummy) ))
                if self.TWIN_SPACE:
                    mat = op.OP1D_tt_fermi(self.mat_creation.T, self.I_elec_sign_R_list, self.I_elec_list, 2*i_state+1, 2*self.n_elec*self.n_spin)
                else:
                    mat = op.OP1D_tt_fermi(self.mat_annihilation_R, self.I_elec_sign_R_list, self.I_elec_list, i_state, self.n_elec*self.n_spin)
                d_creation_R.append( tt.kron(mat, tt.kron(self.II_vib,self.II_dummy) ))
    
                if self.TWIN_SPACE:
                    mat = op.OP1D_tt_fermi(self.mat_annihilation, self.I_elec_sign_L_list, self.I_elec_list, 2*i_state, 2*self.n_elec*self.n_spin)
                else:
                    mat = op.OP1D_tt_fermi(self.mat_annihilation_L, self.I_elec_sign_L_list, self.I_elec_list, i_state, self.n_elec*self.n_spin)
                d_annihilation_L.append( tt.kron(mat, tt.kron(self.II_vib,self.II_dummy) ))
                if self.TWIN_SPACE:
                    mat = op.OP1D_tt_fermi(self.mat_annihilation.T, self.I_elec_sign_R_list, self.I_elec_list, 2*i_state+1, 2*self.n_elec*self.n_spin)
                else:
                    mat = op.OP1D_tt_fermi(self.mat_creation_R, self.I_elec_sign_R_list, self.I_elec_list, i_state, self.n_elec*self.n_spin)
                d_annihilation_R.append( tt.kron(mat, tt.kron(self.II_vib,self.II_dummy) ))
    
                if self.SPIN_ALT_ORDER:
                    i_state = (self.n_elec-1)*self.n_spin + i_spin
                else:
                    i_state = i_spin*self.n_elec + (self.n_elec-1)
                if self.TWIN_SPACE:
                    mat = op.OP1D_tt_fermi(self.mat_creation, self.I_elec_sign_L_list, self.I_elec_list, 2*i_state, 2*self.n_elec*self.n_spin)
                else:
                    mat = op.OP1D_tt_fermi(self.mat_creation_L, self.I_elec_sign_L_list, self.I_elec_list, i_state, self.n_elec*self.n_spin)
                d_creation_L.append( tt.kron(mat, tt.kron(self.II_vib,self.II_dummy) ))
                if self.TWIN_SPACE:
                    mat = op.OP1D_tt_fermi(self.mat_creation.T, self.I_elec_sign_R_list, self.I_elec_list, 2*i_state+1, 2*self.n_elec*self.n_spin)
                else:
                    mat = op.OP1D_tt_fermi(self.mat_annihilation_R, self.I_elec_sign_R_list, self.I_elec_list, i_state, self.n_elec*self.n_spin)
                d_creation_R.append( tt.kron(mat, tt.kron(self.II_vib,self.II_dummy) ))
    
                if self.TWIN_SPACE:
                    mat = op.OP1D_tt_fermi(self.mat_annihilation, self.I_elec_sign_L_list, self.I_elec_list, 2*i_state, 2*self.n_elec*self.n_spin)
                else:
                    mat = op.OP1D_tt_fermi(self.mat_annihilation_L, self.I_elec_sign_L_list, self.I_elec_list, i_state, self.n_elec*self.n_spin)
                d_annihilation_L.append( tt.kron(mat, tt.kron(self.II_vib,self.II_dummy) ))
                if self.TWIN_SPACE:
                    mat = op.OP1D_tt_fermi(self.mat_annihilation.T, self.I_elec_sign_R_list, self.I_elec_list, 2*i_state+1, 2*self.n_elec*self.n_spin)
                else:
                    mat = op.OP1D_tt_fermi(self.mat_creation_R, self.I_elec_sign_R_list, self.I_elec_list, i_state, self.n_elec*self.n_spin)
                d_annihilation_R.append( tt.kron(mat, tt.kron(self.II_vib,self.II_dummy) ))

            self.d_sys_L = []
            self.d_sys_L.append( d_creation_L )
            self.d_sys_L.append( d_annihilation_L )
            self.d_sys_R = []
            self.d_sys_R.append( d_creation_R )
            self.d_sys_R.append( d_annihilation_R )

            # The number of effective environmental electronic DoFs
            self.N_env_imagtime = self.n_spin * self.n_lead * self.n_pole_fermi * 2 * 2
    
            # Identity operators for environmental electronic DoFs
            I_fermibath_list_imagtime=[]
            I_fermibath_sign_list_imagtime=[]
            self.II_fermibath_imagtime=None
            for i_spin in range(self.n_spin):
                for i_lead in range(self.n_lead):
                    for i_pole in range(self.n_pole_fermi):
                        for i_sigma in range(2):
                            for i_sincos in range(2):
                                I_fermibath_list_imagtime.append( tt.eye( [2] ) )
                                I_fermibath_sign_list_imagtime.append( tt.matrix(np.array([[1,0],[0,-1]])) )
                                self.II_fermibath_imagtime = tt.kron( self.II_fermibath_imagtime, tt.eye( [2] ) )

            i_pos = -1
            self.Bp_imagtime = []
            self.Bm_imagtime = []
            for i_spin in range(self.n_spin):
                for i_lead in range(self.n_lead):
                    for i_pole in range(self.n_pole_fermi):
                        for i_sigma in range(2):
                            for i_sincos in range(2):
                                i_pos += 1
                                self.Bp_imagtime.append( 1j * op.OP1D_tt_fermi(self.mat_annihilation, I_fermibath_sign_list_imagtime, I_fermibath_list_imagtime, i_pos, self.N_env_imagtime))
                                self.Bm_imagtime.append( -1j * op.OP1D_tt_fermi(self.mat_creation, I_fermibath_sign_list_imagtime, I_fermibath_list_imagtime, i_pos, self.N_env_imagtime) )

        else:
            d_creation_L = []
            d_creation_R = []
            d_annihilation_L = []
            d_annihilation_R = []
            for i_elec in range(self.n_elec):
                d_creation_L_sub=[]
                d_creation_R_sub=[]
                d_annihilation_L_sub=[]
                d_annihilation_R_sub=[]
                for i_spin in range(self.n_spin):
                    if self.SPIN_ALT_ORDER:
                        i_state = i_elec*self.n_spin + i_spin
                    else:
                        i_state = i_spin*self.n_elec + i_elec
                    if self.TWIN_SPACE:
                        mat = op.OP1D_tt_fermi(self.mat_creation, self.I_elec_sign_L_list, self.I_elec_list, 2*i_state, 2*self.n_elec*self.n_spin)
                    else:
                        mat = op.OP1D_tt_fermi(self.mat_creation_L, self.I_elec_sign_L_list, self.I_elec_list, i_state, self.n_elec*self.n_spin)
                    d_creation_L_sub.append( tt.kron(mat, self.II_vib) )
                    if self.TWIN_SPACE:
                        mat = op.OP1D_tt_fermi(self.mat_creation.T, self.I_elec_sign_R_list, self.I_elec_list, 2*i_state+1, 2*self.n_elec*self.n_spin)
                    else:
                        mat = op.OP1D_tt_fermi(self.mat_annihilation_R, self.I_elec_sign_R_list, self.I_elec_list, i_state, self.n_elec*self.n_spin)
                    d_creation_R_sub.append( tt.kron(mat, self.II_vib) )
    
                    if self.TWIN_SPACE:
                        mat = op.OP1D_tt_fermi(self.mat_annihilation, self.I_elec_sign_L_list, self.I_elec_list, 2*i_state, 2*self.n_elec*self.n_spin)
                    else:
                        mat = op.OP1D_tt_fermi(self.mat_annihilation_L, self.I_elec_sign_L_list, self.I_elec_list, i_state, self.n_elec*self.n_spin)
                    d_annihilation_L_sub.append( tt.kron(mat, self.II_vib) )
                    if self.TWIN_SPACE:
                        mat = op.OP1D_tt_fermi(self.mat_annihilation.T, self.I_elec_sign_R_list, self.I_elec_list, 2*i_state+1, 2*self.n_elec*self.n_spin)
                    else:
                        mat = op.OP1D_tt_fermi(self.mat_creation_R, self.I_elec_sign_R_list, self.I_elec_list, i_state, self.n_elec*self.n_spin)
                    d_annihilation_R_sub.append( tt.kron(mat, self.II_vib) )
                d_creation_L.append(d_creation_L_sub)
                d_creation_R.append(d_creation_R_sub)
                d_annihilation_L.append(d_annihilation_L_sub)
                d_annihilation_R.append(d_annihilation_R_sub)

            self.d_sys_L = []
            self.d_sys_L.append( d_creation_L )
            self.d_sys_L.append( d_annihilation_L )
            self.d_sys_R = []
            self.d_sys_R.append( d_creation_R )
            self.d_sys_R.append( d_annihilation_R )

            # The number of effective environmental electronic DoFs
            self.N_env_imagtime = self.n_spin * self.n_elec * self.n_lead * self.n_pole_fermi * 2 * 2
    
            # Identity operators for environmental electronic DoFs
            I_fermibath_list_imagtime=[]
            I_fermibath_sign_list_imagtime=[]
            self.II_fermibath_imagtime=None
            for i_elec in range(self.n_elec):
                for i_spin in range(self.n_spin):
                    for i_lead in range(self.n_lead):
                        for i_pole in range(self.n_pole_fermi):
                            for i_sigma in range(2):
                                for i_sincos in range(2):
                                    I_fermibath_list_imagtime.append( tt.eye( [2] ) )
                                    I_fermibath_sign_list_imagtime.append( tt.matrix(np.array([[1,0],[0,-1]])) )
                                    self.II_fermibath_imagtime = tt.kron( self.II_fermibath_imagtime, tt.eye( [2] ) )
    
            i_pos = -1
            self.Bp_imagtime = []
            self.Bm_imagtime = []
            for i_elec in range(self.n_elec):
                for i_spin in range(self.n_spin):
                    for i_lead in range(self.n_lead):
                        for i_pole in range(self.n_pole_fermi):
                            for i_sigma in range(2):
                                for i_sincos in range(2):
                                    i_pos += 1
                                    self.Bp_imagtime.append( 1j * op.OP1D_tt_fermi(self.mat_annihilation, I_fermibath_sign_list_imagtime, I_fermibath_list_imagtime, i_pos, self.N_env_imagtime) )
                                    self.Bm_imagtime.append( -1j * op.OP1D_tt_fermi(self.mat_creation, I_fermibath_sign_list_imagtime, I_fermibath_list_imagtime, i_pos, self.N_env_imagtime) )

    def OP_observables(self, EL_STATES_OUTPUT, VIB_STATES_OUTPUT, POP_EL_OUTPUT = None, POP_VIB_OUTPUT = None, N_VIB_OUTPUT = False):

        # Tensors for calculating observables  
        if self.TWIN_SPACE:

            Pop_el = []
            if POP_EL_OUTPUT == 'all':
                for state in itertools.product(np.array([0,1]), repeat = self.n_elec * self.n_spin):
                    OpMat = None  
                    for i_elec in range(self.n_elec * self.n_spin):
                        mat = np.zeros((2,2), dtype=float)
                        i_pop = state[i_elec]
                        mat[i_pop,i_pop] = 1
                        mat = tt.matrix(mat, self.eps)
                        OpMat = tt.kron(OpMat, tt.kron(mat, tt.eye([2])))
                    OpMat = tt.kron(OpMat, tt.kron(self.II_vib, self.II_dummy))
                    OpMat = tt.matvec(OpMat, self.I1_sys)
                    OpMat = OpMat.round(self.eps)
                    Pop_el.append( OpMat )
            elif (POP_EL_OUTPUT != None):
                for state in EL_STATES_OUTPUT:
                    OpMat = None  
                    for i_elec in range(self.n_elec * self.n_spin):
                        mat = np.zeros((2,2), dtype=float)
                        i_pop = int(state[i_elec])
                        mat[i_pop, i_pop] = 1
                        mat = tt.matrix(mat, self.eps)
                        OpMat = tt.kron(OpMat, tt.kron(mat, tt.eye([2])))
                    OpMat = tt.kron(OpMat, tt.kron(self.II_vib,self.II_dummy))
                    OpMat = tt.matvec(OpMat, self.I1_sys)
                    OpMat= OpMat.round(self.eps)
                    Pop_el.append( OpMat )
    
            Pop_vib = []
            if POP_VIB_OUTPUT == 'all':
                for state in itertools.product(np.arange(self.vib_basis[0]), repeat=self.n_vib):
                    OpMat = None  
                    for i_vib in range(self.n_vib):
                        mat = np.zeros((self.vib_basis[i_vib], self.vib_basis[i_vib]), dtype=float)
                        i_pop = state[i_vib]
                        mat[i_pop, i_pop] = 1
                        mat = tt.matrix(mat, self.eps)
                        OpMat = tt.kron(OpMat, tt.kron(mat, tt.eye([self.vib_basis[i_vib]])))
                    OpMat = tt.kron(OpMat,self.II_dummy)
                    OpMat = tt.kron(self.II_elec, OpMat)
                    OpMat = tt.matvec(OpMat, self.I1_sys)
                    OpMat = OpMat.round(self.eps)
                    Pop_vib.append( OpMat )
            elif (POP_VIB_OUTPUT != None):
                for state in VIB_STATES_OUTPUT:
                    OpMat = None  
                    for i_vib in range(self.n_vib):
                        mat = np.zeros((self.vib_basis[i_vib], self.vib_basis[i_vib]), dtype=float)
                        i_pop = int(state[i_vib])
                        mat[i_pop, i_pop] = 1
                        mat = tt.matrix(mat, self.eps)
                        OpMat = tt.kron(OpMat, tt.kron(mat, tt.eye([self.vib_basis[i_vib]])))
                    OpMat = tt.kron(OpMat,self.II_dummy)
                    OpMat = tt.kron(self.II_elec, OpMat)
                    OpMat = tt.matvec(OpMat, self.I1_sys)
                    OpMat = OpMat.round(self.eps)
                    Pop_vib.append( OpMat )
    
            N_vib = []
            if N_VIB_OUTPUT:
                for i_vib in range(self.n_vib):
                    vib_matrix = np.diag( np.arange( self.vib_basis[i_vib] ) )
                    vib_matrix = tt.matrix( vib_matrix, self.eps)
                    OpMat = op.OP1D_tt(vib_matrix, self.I_vib_list, 2*i_vib, 2*self.n_vib)

                    OpMat = tt.kron(self.II_elec, tt.kron(OpMat, self.II_dummy))
                    OpMat = OpMat.round(self.eps)
                    N_vib.append( tt.matvec(OpMat, self.I1_sys) )

            N_site = []
            for i_elec in range(self.n_elec * self.n_spin):
                mat = np.zeros((2,2), dtype=float)
                mat[1,1] = 1
                mat = tt.matrix(mat, self.eps)
                OpMat = op.OP1D_tt(mat, self.I_elec_list, 2*i_elec, 2*self.n_elec * self.n_spin)
                OpMat = tt.kron(OpMat, tt.kron(self.II_vib, self.II_dummy))
                OpMat = OpMat.round(self.eps)
                N_site.append( tt.matvec(OpMat, self.I1_sys) )

    
        return self.I1_sys, Pop_el, Pop_vib, N_vib, N_site

    def OP_hamiltonian_CISS(self, Ei, Uij, Tij, SOC, Tij_Ph, SOC_Ph, vib_freq, vib_elec_coup):

        if self.TWIN_SPACE:

            # System Hamiltonian (Electronic part)
            # H = \sum_{i,s} epsilon_{i,s} d_{i,s}^{dagger}d_{i,s} 
            #   + \sum_{i < j, s, s'} (t_{i,s,j,s'} d_{i,s}^{dagger}d_{j,s'} +h.c.) 
            #   + \sum_{i < j, s, s'} U_{ij} n_{i,s}n_{j,s'}
            H_elec_L = None
            H_elec_R = None
            for i_elec in range(self.n_elec):
                for i_spin in range(self.n_spin):
                    if self.SPIN_ALT_ORDER:
                        i_state = i_elec*self.n_spin + i_spin
                    else:
                        i_state = i_spin*self.n_elec + i_elec
                    mat = op.OP1D_tt(self.mat_occupation, self.I_elec_list, 2*i_state, 2*self.n_elec*self.n_spin)
                    H_elec_L += Ei[i_elec][i_spin]*mat
                    H_elec_L = H_elec_L.round(self.eps)
                    mat = op.OP1D_tt(self.mat_occupation.T, self.I_elec_list, 2*i_state+1, 2*self.n_elec*self.n_spin)
                    H_elec_R += Ei[i_elec][i_spin]*mat
                    H_elec_R = H_elec_R.round(self.eps)
    
            for i_elec in range(self.n_elec):
                for i_spin in range(self.n_spin):
                    if self.SPIN_ALT_ORDER:
                        i_state = i_elec*self.n_spin + i_spin
                    else:
                        i_state = i_spin*self.n_elec + i_elec
                    mat_1_L = op.OP1D_tt(self.mat_occupation, self.I_elec_list, 2*i_state, 2*self.n_elec*self.n_spin)
                    mat_1_R = op.OP1D_tt(self.mat_occupation.T, self.I_elec_list, 2*i_state+1, 2*self.n_elec*self.n_spin)
                    for j_elec in range(self.n_elec):
                        for j_spin in range(self.n_spin):
                            if self.SPIN_ALT_ORDER:
                                j_state = j_elec*self.n_spin + j_spin
                            else:
                                j_state = j_spin*self.n_elec + j_elec
                            if (j_state > i_state):
                                mat_2_L = op.OP1D_tt(self.mat_occupation, self.I_elec_list, 2*j_state, 2*self.n_elec*self.n_spin)
                                mat_2_R = op.OP1D_tt(self.mat_occupation.T, self.I_elec_list, 2*j_state+1, 2*self.n_elec*self.n_spin)
                                H_elec_L += Uij[i_elec, i_spin, j_elec, j_spin] * mat_1_L * mat_2_L
                                H_elec_L = H_elec_L.round(self.eps)
                                H_elec_R += Uij[i_elec, i_spin, j_elec, j_spin] * mat_1_R * mat_2_R 
                                H_elec_R = H_elec_R.round(self.eps)
    
            for i_elec in range(self.n_elec):
                for i_spin in range(self.n_spin):
                    if self.SPIN_ALT_ORDER:
                        i_state = i_elec*self.n_spin + i_spin
                    else:
                        i_state = i_spin*self.n_elec + i_elec
                    mat_1_L = op.OP1D_tt_fermi(self.mat_creation, self.I_elec_sign_L_list, self.I_elec_list, 2*i_state, 2*self.n_elec*self.n_spin)
                    mat_1_R = op.OP1D_tt_fermi(self.mat_creation.T, self.I_elec_sign_R_list, self.I_elec_list, 2*i_state+1, 2*self.n_elec*self.n_spin)
                    for j_elec in range(self.n_elec):
                        for j_spin in range(self.n_spin):
                            if self.SPIN_ALT_ORDER:
                                j_state = j_elec*self.n_spin + j_spin
                            else:
                                j_state = j_spin*self.n_elec + j_elec
                            if (j_state > i_state):
                                mat_2_L = op.OP1D_tt_fermi(self.mat_annihilation, self.I_elec_sign_L_list, self.I_elec_list, 2*j_state, 2*self.n_elec*self.n_spin)
                                mat_2_R = op.OP1D_tt_fermi(self.mat_annihilation.T, self.I_elec_sign_R_list, self.I_elec_list, 2*j_state+1, 2*self.n_elec*self.n_spin)
                                mat = mat_1_L * mat_2_L
                                H_elec_L += Tij[i_elec, i_spin, j_elec, j_spin]*(mat + mat.T)
                                H_elec_L = H_elec_L.round(self.eps)
                                mat = mat_1_R * mat_2_R
                                H_elec_R -= Tij[i_elec, i_spin, j_elec, j_spin]*(mat + mat.T)
                                H_elec_R = H_elec_R.round(self.eps)
            
            # Fransson's next-nearest-neighbor spin-orbit coupling
            for i_elec in range(self.n_elec-2):
                for i_spin in range(self.n_spin):
                    if self.SPIN_ALT_ORDER:
                        i_state = i_elec * self.n_spin + i_spin
                    else:
                        i_state = i_spin * self.n_elec + i_elec
                    mat_1_L = op.OP1D_tt_fermi(self.mat_creation, self.I_elec_sign_L_list, self.I_elec_list, 2*i_state, 2*self.n_elec*self.n_spin)
                    mat_1_R = op.OP1D_tt_fermi(self.mat_creation.T, self.I_elec_sign_R_list, self.I_elec_list, 2*i_state+1, 2*self.n_elec*self.n_spin)
                    for j_spin in range(self.n_spin):
                        j_elec = i_elec + 2
                        if self.SPIN_ALT_ORDER:
                            j_state = j_elec * self.n_spin + j_spin
                        else:
                            j_state = j_spin * self.n_elec + j_elec
                        
                        mat_2_L = op.OP1D_tt_fermi(self.mat_annihilation, self.I_elec_sign_L_list, self.I_elec_list, 2*j_state, 2*self.n_elec*self.n_spin)
                        mat_2_R = op.OP1D_tt_fermi(self.mat_annihilation.T, self.I_elec_sign_R_list, self.I_elec_list, 2*j_state+1, 2*self.n_elec*self.n_spin)
                    
                        mat = mat_1_L * mat_2_L
                        H_elec_L += SOC[i_elec, i_spin, j_spin]*(mat + mat.T)
                        H_elec_L = H_elec_L.round(self.eps)
                        mat = mat_1_R * mat_2_R
                        H_elec_R -= SOC[i_elec, i_spin, j_spin]*(mat + mat.T)
                        H_elec_R = H_elec_R.round(self.eps)
             
            self.Liouville_L = tt.kron(H_elec_L, tt.kron(self.II_vib,self.II_dummy))
            self.Liouville_R = tt.kron(H_elec_R, tt.kron(self.II_vib,self.II_dummy))
    
            # System Hamiltonian (Vibrational  part)
            # H = \sum_{i} omega_i a_i^{dagger}a_i 
            H_vib_L = None
            H_vib_R = None
            for i_vib in range(self.n_vib):
                vib_matrix = np.diag( vib_freq[i_vib] * np.arange( self.vib_basis[i_vib] ) )
                vib_matrix = tt.matrix( vib_matrix, self.eps)
                H_vib_L += op.OP1D_tt(vib_matrix, self.I_vib_list, 2*i_vib, 2*self.n_vib)
                H_vib_L = H_vib_L.round(self.eps)
                H_vib_R += op.OP1D_tt(vib_matrix.T, self.I_vib_list, 2*i_vib+1, 2*self.n_vib)
                H_vib_R = H_vib_R.round(self.eps)
            self.Liouville_L += tt.kron(self.II_elec, tt.kron(H_vib_L,self.II_dummy))
            self.Liouville_R += tt.kron(self.II_elec, tt.kron(H_vib_R,self.II_dummy))
    
            # System Hamiltonian (Vibronic coupling)
            Coup_elec_L = []
            Coup_elec_R = []
            for i_elec in range(self.n_elec):
                for i_spin in range(self.n_spin):
                    if self.SPIN_ALT_ORDER:
                        i_state = i_elec*self.n_spin + i_spin
                    else:
                        i_state = i_spin*self.n_elec + i_elec
                    Coup_elec_L.append( op.OP1D_tt(self.mat_occupation, self.I_elec_list, 2*i_state, 2*self.n_elec*self.n_spin) )
                    Coup_elec_R.append( op.OP1D_tt(self.mat_occupation.T, self.I_elec_list, 2*i_state+1, 2*self.n_elec*self.n_spin) )
            
            
            self.Coup_vib_L = []
            self.Coup_vib_R = []
            for i_vib in range(self.n_vib):
                mat_vib = np.diag(np.sqrt(np.arange(1,self.vib_basis[i_vib])),k=-1)
                mat_vib += mat_vib.T
                mat_vib = tt.matrix( mat_vib, self.eps )
                self.Coup_vib_L.append( op.OP1D_tt(mat_vib, self.I_vib_list, 2*i_vib, 2*self.n_vib) )
                self.Coup_vib_R.append( op.OP1D_tt(mat_vib.T, self.I_vib_list, 2*i_vib+1, 2*self.n_vib) )
    
            Coup_elecvib_L = None
            Coup_elecvib_R = None
            for i_elec in range(self.n_elec):
                for i_spin in range(self.n_spin):
                    if self.SPIN_ALT_ORDER:
                        i_state = i_elec*self.n_spin + i_spin
                    else:
                        i_state = i_spin*self.n_elec + i_elec
                    for i_vib in range(self.n_vib):
                        Coup_elecvib_L += vib_elec_coup[i_elec][i_spin][i_vib] * tt.kron(Coup_elec_L[i_state], tt.kron(self.Coup_vib_L[i_vib],self.II_dummy))
                        Coup_elecvib_L = Coup_elecvib_L.round(self.eps)
                        Coup_elecvib_R += vib_elec_coup[i_elec][i_spin][i_vib] * tt.kron(Coup_elec_R[i_state], tt.kron(self.Coup_vib_R[i_vib],self.II_dummy))
                        Coup_elecvib_R = Coup_elecvib_R.round(self.eps)
            self.Liouville_L += Coup_elecvib_L
            self.Liouville_R += Coup_elecvib_R
    
            # Fransson's nearest-neighbor electron-phonon coupling (spin-independent)
            Coup_elecvib_nn_L = None
            Coup_elecvib_nn_R = None
            for i_elec in range(self.n_elec-1):
                for i_spin in range(self.n_spin):
                    if self.SPIN_ALT_ORDER:
                        i_state = i_elec*self.n_spin + i_spin
                    else:
                        i_state = i_spin*self.n_elec + i_elec
                    mat_1_L = op.OP1D_tt_fermi(self.mat_creation, self.I_elec_sign_L_list, self.I_elec_list, 2*i_state, 2*self.n_elec*self.n_spin)
                    mat_1_R = op.OP1D_tt_fermi(self.mat_creation.T, self.I_elec_sign_R_list, self.I_elec_list, 2*i_state+1, 2*self.n_elec*self.n_spin)     
    
                    for j_spin in range(self.n_spin):
                        j_elec = i_elec + 1
            
                        if self.SPIN_ALT_ORDER:
                            j_state = j_elec*self.n_spin + j_spin
                        else:
                            j_state = j_spin*self.n_elec + j_elec
                        
                            
                        mat_2_L = op.OP1D_tt_fermi(self.mat_annihilation, self.I_elec_sign_L_list, self.I_elec_list, 2*j_state, 2*self.n_elec*self.n_spin)
                        mat_2_R = op.OP1D_tt_fermi(self.mat_annihilation.T, self.I_elec_sign_R_list, self.I_elec_list, 2*j_state+1, 2*self.n_elec*self.n_spin)
            
                        mat = mat_1_L * mat_2_L
                        Coup_elecvib_nn_L += Tij_Ph[i_elec, i_spin, j_spin] * tt.kron((mat+mat.T), tt.kron(self.Coup_vib_L[i_elec],self.II_dummy))
                        Coup_elecvib_nn_L = Coup_elecvib_nn_L.round(self.eps)
                        mat = mat_1_R * mat_2_R
                        # (?) either + or minus
                        Coup_elecvib_nn_R += -Tij_Ph[i_elec, i_spin, j_spin] * tt.kron((mat+mat.T), tt.kron(self.Coup_vib_R[i_elec],self.II_dummy))
                        Coup_elecvib_nn_R = Coup_elecvib_nn_R.round(self.eps)
                        
            self.Liouville_L += Coup_elecvib_nn_L
            self.Liouville_R += Coup_elecvib_nn_R
    
    
            # Fransson's next-nearest-neighbor spin-dependent electron-phonon coupling
            Coup_elecvib_soc_L = None
            Coup_elecvib_soc_R = None
            for i_elec in range(self.n_elec-2):
                for i_spin in range(self.n_spin):
                    if self.SPIN_ALT_ORDER:
                        i_state = i_elec*self.n_spin + i_spin
                    else:
                        i_state = i_spin*self.n_elec + i_elec
                    mat_1_L = op.OP1D_tt_fermi(self.mat_creation, self.I_elec_sign_L_list, self.I_elec_list, 2*i_state, 2*self.n_elec*self.n_spin)
                    mat_1_R = op.OP1D_tt_fermi(self.mat_creation.T, self.I_elec_sign_R_list, self.I_elec_list, 2*i_state+1, 2*self.n_elec*self.n_spin)     
                    
                    for j_spin in range(self.n_spin):
                        j_elec = i_elec + 2
                        if self.SPIN_ALT_ORDER:
                            j_state = j_elec*self.n_spin + j_spin
                        else:
                            j_state = j_spin*self.n_elec + j_elec
                        
                            
                        mat_2_L = op.OP1D_tt_fermi(self.mat_annihilation, self.I_elec_sign_L_list, self.I_elec_list, 2*j_state, 2*self.n_elec*self.n_spin)
                        mat_2_R = op.OP1D_tt_fermi(self.mat_annihilation.T, self.I_elec_sign_R_list, self.I_elec_list, 2*j_state+1, 2*self.n_elec*self.n_spin)

                        mat = mat_1_L * mat_2_L
                        Coup_elecvib_soc_L += SOC_Ph[i_elec, i_spin, j_spin] * tt.kron((mat+mat.T), tt.kron(self.Coup_vib_L[i_elec],self.II_dummy))
                        Coup_elecvib_soc_L = Coup_elecvib_soc_L.round(self.eps)
                        
                        mat = mat_1_R * mat_2_R
                        Coup_elecvib_soc_R += -SOC_Ph[i_elec, i_spin, j_spin] * tt.kron((mat+mat.T), tt.kron(self.Coup_vib_R[i_elec],self.II_dummy))
    
                        Coup_elecvib_soc_R = Coup_elecvib_soc_R.round(self.eps)
            self.Liouville_L += Coup_elecvib_soc_L
            self.Liouville_R += Coup_elecvib_soc_R
            self.Liouville_L = self.Liouville_L.round(self.eps)
            self.Liouville_R = self.Liouville_R.round(self.eps)

        else:

            # System Hamiltonian (Electronic part)
            # H = \sum_{i,s} epsilon_{i,s} d_{i,s}^{dagger}d_{i,s} 
            #   + \sum_{i < j, s, s'} (t_{i,s,j,s'} d_{i,s}^{dagger}d_{j,s'} +h.c.) 
            #   + \sum_{i < j, s, s'} U_{ij} n_{i,s}n_{j,s'}
            H_elec_L = None
            H_elec_R = None
            for i_elec in range(self.n_elec):
                for i_spin in range(self.n_spin):
                    if self.SPIN_ALT_ORDER:
                        i_state = i_elec*self.n_spin + i_spin
                    else:
                        i_state = i_spin*self.n_elec + i_elec
                    mat = op.OP1D_tt(self.mat_occupation_L, self.I_elec_list, i_state, self.n_elec*self.n_spin)
                    H_elec_L += Ei[i_elec][i_spin]*mat
                    H_elec_L = H_elec_L.round(self.eps)
                    mat = op.OP1D_tt(self.mat_occupation_R, self.I_elec_list, i_state, self.n_elec*self.n_spin)
                    H_elec_R += Ei[i_elec][i_spin]*mat
                    H_elec_R = H_elec_R.round(self.eps)
    
            for i_elec in range(self.n_elec):
                for i_spin in range(self.n_spin):
                    if self.SPIN_ALT_ORDER:
                        i_state = i_elec*self.n_spin + i_spin
                    else:
                        i_state = i_spin*self.n_elec + i_elec
                    mat_1_L = op.OP1D_tt(self.mat_occupation_L, self.I_elec_list, i_state, self.n_elec*self.n_spin)
                    mat_1_R = op.OP1D_tt(self.mat_occupation_R, self.I_elec_list, i_state, self.n_elec*self.n_spin)
                    for j_elec in range(self.n_elec):
                        for j_spin in range(self.n_spin):
                            if self.SPIN_ALT_ORDER:
                                j_state = j_elec*self.n_spin + j_spin
                            else:
                                j_state = j_spin*self.n_elec + j_elec
                            if (j_state > i_state):
                                mat_2_L = op.OP1D_tt(self.mat_occupation_L, self.I_elec_list, j_state, self.n_elec*self.n_spin)
                                mat_2_R = op.OP1D_tt(self.mat_occupation_R, self.I_elec_list, j_state, self.n_elec*self.n_spin)
                                H_elec_L += Uij[i_elec, i_spin, j_elec, j_spin] * mat_1_L * mat_2_L
                                H_elec_L = H_elec_L.round(self.eps)
                                H_elec_R += Uij[i_elec, i_spin, j_elec, j_spin] * mat_1_R * mat_2_R 
                                H_elec_R = H_elec_R.round(self.eps)
    
            for i_elec in range(self.n_elec):
                for i_spin in range(self.n_spin):
                    if self.SPIN_ALT_ORDER:
                        i_state = i_elec*self.n_spin + i_spin
                    else:
                        i_state = i_spin*self.n_elec + i_elec
                    mat_1_L = op.OP1D_tt_fermi(self.mat_creation_L, self.I_elec_sign_L_list, self.I_elec_list, i_state, self.n_elec*self.n_spin)
                    mat_1_R = op.OP1D_tt_fermi(self.mat_annihilation_R, self.I_elec_sign_R_list, self.I_elec_list, i_state, self.n_elec*self.n_spin)
                    for j_elec in range(self.n_elec):
                        for j_spin in range(self.n_spin):
                            if self.SPIN_ALT_ORDER:
                                j_state = j_elec*self.n_spin + j_spin
                            else:
                                j_state = j_spin*self.n_elec + j_elec
                            if (j_state > i_state):
                                mat_2_L = op.OP1D_tt_fermi(self.mat_annihilation_L, self.I_elec_sign_L_list, self.I_elec_list, j_state, self.n_elec*self.n_spin)
                                mat_2_R = op.OP1D_tt_fermi(self.mat_creation_R, self.I_elec_sign_R_list, self.I_elec_list, j_state, self.n_elec*self.n_spin)
                                mat = mat_1_L * mat_2_L
                                H_elec_L += Tij[i_elec, i_spin, j_elec, j_spin]*(mat + mat.T)
                                H_elec_L = H_elec_L.round(self.eps)
                                mat = mat_1_R * mat_2_R
                                H_elec_R -= Tij[i_elec, i_spin, j_elec, j_spin]*(mat + mat.T)
                                H_elec_R = H_elec_R.round(self.eps)
            
            # Fransson's next-nearest-neighbor spin-orbit coupling
            for i_elec in range(self.n_elec-2):
                for i_spin in range(self.n_spin):
                    if self.SPIN_ALT_ORDER:
                        i_state = i_elec * self.n_spin + i_spin
                    else:
                        i_state = i_spin * self.n_elec + i_elec
                    mat_1_L = op.OP1D_tt_fermi(self.mat_creation_L, self.I_elec_sign_L_list, self.I_elec_list, i_state, self.n_elec*self.n_spin)
                    mat_1_R = op.OP1D_tt_fermi(self.mat_annihilation_R, self.I_elec_sign_R_list, self.I_elec_list, i_state, self.n_elec*self.n_spin)
                    for j_spin in range(self.n_spin):
                        j_elec = i_elec + 2
                        if self.SPIN_ALT_ORDER:
                            j_state = j_elec * self.n_spin + j_spin
                        else:
                            j_state = j_spin * self.n_elec + j_elec
                        
                        mat_2_L = op.OP1D_tt_fermi(self.mat_annihilation_L, self.I_elec_sign_L_list, self.I_elec_list, j_state, self.n_elec*self.n_spin)
                        mat_2_R = op.OP1D_tt_fermi(self.mat_creation_R, self.I_elec_sign_R_list, self.I_elec_list, j_state, self.n_elec*self.n_spin)
                    
                        mat = mat_1_L * mat_2_L
                        H_elec_L += SOC[i_elec, i_spin, j_spin]*(mat + mat.T)
                        H_elec_L = H_elec_L.round(self.eps)
                        mat = mat_1_R * mat_2_R
                        H_elec_R -= SOC[i_elec, i_spin, j_spin]*(mat + mat.T)
                        H_elec_R = H_elec_R.round(self.eps)
             
            self.Liouville_L = tt.kron(H_elec_L, self.II_vib)
            self.Liouville_R = tt.kron(H_elec_R, self.II_vib)
    
            # System Hamiltonian (Vibrational  part)
            # H = \sum_{i} omega_i a_i^{dagger}a_i 
            H_vib_L = None
            H_vib_R = None
            for i_vib in range(self.n_vib):
                vib_matrix = np.diag( vib_freq[i_vib] * np.arange( self.vib_basis[i_vib] ) )
                vib_matrix_L = np.kron( vib_matrix, np.eye(self.vib_basis[i_vib]) )
                vib_matrix_L = tt.matrix( vib_matrix_L, self.eps)
                vib_matrix_R = np.kron( np.eye(self.vib_basis[i_vib]), vib_matrix.T )
                vib_matrix_R = tt.matrix( vib_matrix_R, self.eps)
                H_vib_L += op.OP1D_tt(vib_matrix_L, self.I_vib_list, i_vib, self.n_vib)
                H_vib_L = H_vib_L.round(self.eps)
                H_vib_R += op.OP1D_tt(vib_matrix_R, self.I_vib_list, i_vib, self.n_vib)
                H_vib_R = H_vib_R.round(self.eps)
            self.Liouville_L += tt.kron(self.II_elec, H_vib_L)
            self.Liouville_R += tt.kron(self.II_elec, H_vib_R)
    
            # System Hamiltonian (Vibronic coupling)
            Coup_elec_L = []
            Coup_elec_R = []
            for i_elec in range(self.n_elec):
                for i_spin in range(self.n_spin):
                    if self.SPIN_ALT_ORDER:
                        i_state = i_elec*self.n_spin + i_spin
                    else:
                        i_state = i_spin*self.n_elec + i_elec
                    Coup_elec_L.append( op.OP1D_tt(self.mat_occupation_L, self.I_elec_list, i_state, self.n_elec*self.n_spin) )
                    Coup_elec_R.append( op.OP1D_tt(self.mat_occupation_R, self.I_elec_list, i_state, self.n_elec*self.n_spin) )
            
            
            self.Coup_vib_L = []
            self.Coup_vib_R = []
            for i_vib in range(self.n_vib):
                mat_vib = np.diag(np.sqrt(np.arange(1,self.vib_basis[i_vib]))/np.sqrt(2.),k=-1)
                mat_vib += mat_vib.T
                mat_vib_L = np.kron( mat_vib, np.eye(self.vib_basis[i_vib]) )
                mat_vib_L = tt.matrix( mat_vib_L, self.eps )
                mat_vib_R = np.kron( np.eye(self.vib_basis[i_vib]), mat_vib.T )
                mat_vib_R = tt.matrix( mat_vib_R, self.eps )
                self.Coup_vib_L.append( op.OP1D_tt(mat_vib_L, self.I_vib_list, i_vib, self.n_vib) )
                self.Coup_vib_R.append( op.OP1D_tt(mat_vib_R, self.I_vib_list, i_vib, self.n_vib) )
    
            Coup_elecvib_L = None
            Coup_elecvib_R = None
            for i_elec in range(self.n_elec):
                for i_spin in range(self.n_spin):
                    if self.SPIN_ALT_ORDER:
                        i_state = i_elec*self.n_spin + i_spin
                    else:
                        i_state = i_spin*self.n_elec + i_elec
                    for i_vib in range(self.n_vib):
                        Coup_elecvib_L += vib_elec_coup[i_elec][i_spin][i_vib] * tt.kron(Coup_elec_L[i_state], self.Coup_vib_L[i_vib])
                        Coup_elecvib_L = Coup_elecvib_L.round(self.eps)
                        Coup_elecvib_R += vib_elec_coup[i_elec][i_spin][i_vib] * tt.kron(Coup_elec_R[i_state], self.Coup_vib_R[i_vib])
                        Coup_elecvib_R = Coup_elecvib_R.round(self.eps)
            self.Liouville_L += Coup_elecvib_L
            self.Liouville_R += Coup_elecvib_R
    
            # Fransson's nearest-neighbor electron-phonon coupling (spin-independent)
            Coup_elecvib_nn_L = None
            Coup_elecvib_nn_R = None
            for i_elec in range(self.n_elec-1):
                for i_spin in range(self.n_spin):
                    if self.SPIN_ALT_ORDER:
                        i_state = i_elec*self.n_spin + i_spin
                    else:
                        i_state = i_spin*self.n_elec + i_elec
                    mat_1_L = op.OP1D_tt_fermi(self.mat_creation_L, self.I_elec_sign_L_list, self.I_elec_list, i_state, self.n_elec*self.n_spin)
                    mat_1_R = op.OP1D_tt_fermi(self.mat_annihilation_R, self.I_elec_sign_R_list, self.I_elec_list, i_state, self.n_elec*self.n_spin)     
    
                    for j_spin in range(self.n_spin):
                        j_elec = i_elec + 1
            
                        if self.SPIN_ALT_ORDER:
                            j_state = j_elec*self.n_spin + j_spin
                        else:
                            j_state = j_spin*self.n_elec + j_elec
                        
                            
                        mat_2_L = op.OP1D_tt_fermi(self.mat_annihilation_L, self.I_elec_sign_L_list, self.I_elec_list, j_state, self.n_elec*self.n_spin)
                        mat_2_R = op.OP1D_tt_fermi(self.mat_creation_R, self.I_elec_sign_R_list, self.I_elec_list, j_state, self.n_elec*self.n_spin)
            
                        mat = mat_1_L * mat_2_L
                        Coup_elecvib_nn_L += Tij_Ph[i_elec, i_spin, j_spin] * tt.kron((mat+mat.T), self.Coup_vib_L[i_elec])
                        Coup_elecvib_nn_L = Coup_elecvib_nn_L.round(self.eps)
                        mat = mat_1_R * mat_2_R
                        Coup_elecvib_nn_R = -Tij_Ph[i_elec, i_spin, j_spin] * tt.kron((mat+mat.T), self.Coup_vib_R[i_elec])
                        Coup_elecvib_nn_R = Coup_elecvib_nn_R.round(self.eps)
                        
            self.Liouville_L += Coup_elecvib_nn_L
            self.Liouville_R += Coup_elecvib_nn_R
    
    
            # Fransson's next-nearest-neighbor spin-dependent electron-phonon coupling
            Coup_elecvib_soc_L = None
            Coup_elecvib_soc_R = None
            for i_elec in range(self.n_elec-2):
                for i_spin in range(self.n_spin):
                    if self.SPIN_ALT_ORDER:
                        i_state = i_elec*self.n_spin + i_spin
                    else:
                        i_state = i_spin*self.n_elec + i_elec
                    mat_1_L = op.OP1D_tt_fermi(self.mat_creation_L, self.I_elec_sign_L_list, self.I_elec_list, i_state, self.n_elec*self.n_spin)
                    mat_1_R = op.OP1D_tt_fermi(self.mat_annihilation_R, self.I_elec_sign_R_list, self.I_elec_list, i_state, self.n_elec*self.n_spin)     
                    
                    for j_spin in range(self.n_spin):
                        j_elec = i_elec + 2
                        if self.SPIN_ALT_ORDER:
                            j_state = j_elec*self.n_spin + j_spin
                        else:
                            j_state = j_spin*self.n_elec + j_elec
                        
                            
                        mat_2_L = op.OP1D_tt_fermi(self.mat_annihilation_L, self.I_elec_sign_L_list, self.I_elec_list, j_state, self.n_elec*self.n_spin)
                        mat_2_R = op.OP1D_tt_fermi(self.mat_creation_R, self.I_elec_sign_R_list, self.I_elec_list, j_state, self.n_elec*self.n_spin)

                        mat = mat_1_L * mat_2_L
                        Coup_elecvib_soc_L += SOC_Ph[i_elec, i_spin, j_spin] * tt.kron((mat+mat.T), self.Coup_vib_L[i_elec])
                        Coup_elecvib_soc_L = Coup_elecvib_soc_L.round(self.eps)
                        
                        mat = mat_1_R * mat_2_R
                        Coup_elecvib_soc_R += -SOC_Ph[i_elec, i_spin, j_spin] * tt.kron((mat+mat.T), self.Coup_vib_R[i_elec])
    
                        Coup_elecvib_soc_R = Coup_elecvib_soc_R.round(self.eps)
            self.Liouville_L += Coup_elecvib_soc_L
            self.Liouville_R += Coup_elecvib_soc_R
            self.Liouville_L = self.Liouville_L.round(self.eps)
            self.Liouville_R = self.Liouville_R.round(self.eps)

    def HQME_fermi_env_realtime(self, coefficients_fermi, exponents_fermi, elec_lead_couplings, WBL = False):

        if self.HEAD_TAIL_COUP:

            # The number of effective environmental electronic DoFs
            N_env = self.n_spin * self.n_lead * self.n_pole_fermi * 2 + self.n_dummy
    
            # Identity operators for environmental electronic DoFs
            I_fermibath_list=[]
            I_fermibath_sign_list=[]
            II_fermibath=None
            II_fermibath_sign=None
            for i in range(self.n_dummy):
                I_fermibath_list.append(tt.eye([2]))            
                I_fermibath_sign_list.append(tt.eye([2]))
                II_fermibath = tt.kron(II_fermibath, tt.eye([2]))
                II_fermibath_sign = tt.kron(II_fermibath_sign, tt.eye([2]))
            

            for i_spin in range(self.n_spin):
                for i_lead in range(self.n_lead):
                    for i_pole in range(self.n_pole_fermi):
                        for i_sigma in range(2):
                            I_fermibath_list.append( tt.eye( [2] ) )
                            I_fermibath_sign_list.append( tt.matrix(np.array([[1,0],[0,-1]])) )
                            II_fermibath = tt.kron( II_fermibath, tt.eye( [2] ) )
                            II_fermibath_sign = tt.kron( II_fermibath_sign, tt.matrix(np.array([[1,0],[0,-1]])) )
            
            
            #Generate the vector for the bath indices |0, 0, ..., 0> (with N_bath*N_pole zeros)
            v0 = None
            for i in range(N_env):
                vec_temp = np.zeros((2), dtype=complex) 
                vec_temp[0] = 1.+0.j
                vec_temp = tt.vector( vec_temp )
                v0 = tt.kron(v0, vec_temp)

            I1_dummy_R = None
            for i in range(self.n_dummy):
                vec_temp = np.zeros((2), dtype=complex)
                vec_temp[0] = 1.+0.j
                vec_temp = tt.vector( vec_temp )
                I1_dummy_R = tt.kron(I1_dummy_R, vec_temp)

            Ix0=tt.mkron(v0, self.I1_elec, self.I1_vib, self.I1_dummy)
    
            #Generate the vectors for the bath indices |1, 0, ..., 0>, |0, 1, ..., 0>,... 
            v1 = [None]*N_env
            vec_temp1 = tt.vector( np.array([1,0],dtype=complex) )
            vec_temp2 = tt.vector( np.array([0,1],dtype=complex) )
            for i in range(N_env):
                for j in range(N_env):
                    if (j == i):
                        v1[i] = tt.kron(v1[i], vec_temp2)
                    else:
                        v1[i] = tt.kron(v1[i], vec_temp1)
            

            XOP = None
            XOP += -1j*tt.kron(II_fermibath, (self.Liouville_L-self.Liouville_R))
            XOP = XOP.round(self.eps)
            
            # Super-Hamiltonian for system-bath interaction in the HEOM method
            i_pos = -1 + self.n_dummy
            Gop = None
            Sop = None
            Rop = None
            Current = [None] * self.n_spin * self.n_lead
            for i_spin in range(self.n_spin):
                for i_lead in range(self.n_lead):
                    i_state = i_spin*self.n_lead+i_lead
                    if WBL:
                        W1 = tt.kron(II_fermibath, self.d_sys_L[1][i_state]) - tt.kron(II_fermibath_sign, self.d_sys_R[1][i_state])
                        W2 = tt.kron(II_fermibath, self.d_sys_L[0][i_state]) - tt.kron(II_fermibath_sign, self.d_sys_R[0][i_state])
                        XOP += -elec_lead_couplings[i_spin][i_lead]**2/4.* ( W1*W2 + W2*W1)
                        XOP = XOP.round(self.eps)
    
                        OpMat = tt.kron(II_fermibath, self.d_sys_L[0][i_state] * self.d_sys_L[1][i_state])
                        Current[i_state] += -elec_lead_couplings[i_spin][i_lead]**2/2.*tt.matvec(OpMat, Ix0)
                        OpMat = tt.kron(II_fermibath, self.d_sys_L[1][i_state] * self.d_sys_L[0][i_state])
                        Current[i_state] += elec_lead_couplings[i_spin][i_lead]**2/2.*tt.matvec(OpMat, Ix0)
                    for i_pole in range(self.n_pole_fermi):
                        for i_sigma in range(2):
                            i_pos += 1
    
                            BB = op.OP1D_tt(self.mat_occupation, I_fermibath_list, i_pos, N_env)
                            Gop += -exponents_fermi[i_lead][i_spin][i_pole][i_sigma] * tt.kron(BB, self.II_sys)
    
                            Bp = op.OP1D_tt_fermi(self.mat_annihilation, I_fermibath_sign_list, I_fermibath_list, i_pos, N_env)
                            Sop += -1j * elec_lead_couplings[i_spin][i_lead] * tt.kron(Bp, self.d_sys_L[1-i_sigma][i_state])

                            Bp = op.OP1D_tt_fermi(self.mat_annihilation, I_fermibath_list, I_fermibath_sign_list, i_pos, N_env)
                            Sop += 1j * elec_lead_couplings[i_spin][i_lead] * tt.kron(Bp, self.d_sys_R[1-i_sigma][i_state])
                            
                            Bm = op.OP1D_tt_fermi(self.mat_creation, I_fermibath_sign_list, I_fermibath_list, i_pos, N_env)
                            Rop += -1j * elec_lead_couplings[i_spin][i_lead] * coefficients_fermi[i_lead][i_spin][i_pole] * tt.kron(Bm, self.d_sys_L[i_sigma][i_state])

                            Bm = op.OP1D_tt_fermi(self.mat_creation, I_fermibath_list, I_fermibath_sign_list, i_pos, N_env)
                            Rop += 1j * elec_lead_couplings[i_spin][i_lead] * np.conjugate( coefficients_fermi[i_lead][i_spin][i_pole] ) * tt.kron(Bm, self.d_sys_R[i_sigma][i_state])
                            OpMat = tt.kron(II_fermibath, self.d_sys_L[i_sigma][i_state]) 
                            Ix=tt.mkron(v1[i_pos], self.I1_elec, self.I1_vib, self.I1_dummy)
                            Current[i_state] += -1j* (-1)**i_sigma * elec_lead_couplings[i_spin][i_lead] * tt.matvec(OpMat, Ix)
    
                            Gop = Gop.round(self.eps)
                            Sop = Sop.round(self.eps)
                            Rop = Rop.round(self.eps)
                    Current[i_state] = Current[i_state].round(self.eps)
    
        else:

            # The number of effective environmental electronic DoFs
            N_env = self.n_spin * self.n_elec * self.n_lead * self.n_pole_fermi * 2
    
            # Identity operators for environmental electronic DoFs
            I_fermibath_list=[]
            I_fermibath_sign_list=[]
            II_fermibath=None
            II_fermibath_sign=None
            for i_elec in range(self.n_elec):
                for i_spin in range(self.n_spin):
                    for i_lead in range(self.n_lead):
                        for i_pole in range(self.n_pole_fermi):
                            for i_sigma in range(2):
                                I_fermibath_list.append( tt.eye( [2] ) )
                                I_fermibath_sign_list.append( tt.matrix(np.array([[1,0],[0,-1]])) )
                                II_fermibath = tt.kron( II_fermibath, tt.eye( [2] ) )
                                II_fermibath_sign = tt.kron( II_fermibath_sign, tt.matrix(np.array([[1,0],[0,-1]])) )
    
            #Generate the vector for the bath indices |0, 0, ..., 0> (with N_bath*N_pole zeros)
            v0 = None
            for i in range(N_env):
                vec_temp = np.zeros((2), dtype=complex) 
                vec_temp[0] = 1.+0.j
                vec_temp = tt.vector( vec_temp )
                v0 = tt.kron(v0, vec_temp)
            Ix0=tt.mkron(v0, self.I1_elec, self.I1_vib)
    
            #Generate the vectors for the bath indices |1, 0, ..., 0>, |0, 1, ..., 0>,... 
            v1 = [None]*N_env
            vec_temp1 = tt.vector( np.array([1,0],dtype=complex) )
            vec_temp2 = tt.vector( np.array([0,1],dtype=complex) )
            for i in range(N_env):
                for j in range(N_env):
                    if (j == i):
                        v1[i] = tt.kron(v1[i], vec_temp2)
                    else:
                        v1[i] = tt.kron(v1[i], vec_temp1)

            XOP = None
            XOP += -1j*tt.kron(II_fermibath, (self.Liouville_L-self.Liouville_R))
            XOP = XOP.round(self.eps)
    
            i_pos = -1
            Gop = None
            Sop = None
            Rop = None
            Current = [None] * self.n_spin * self.n_lead
            for i_elec in range(self.n_elec):
                for i_spin in range(self.n_spin):
                    for i_lead in range(self.n_lead):
                        i_cstate = i_spin*self.n_lead+i_lead
                        if WBL:
                            W1 = tt.kron(II_fermibath, self.d_sys_L[1][i_elec][i_spin]) - tt.kron(II_fermibath_sign, self.d_sys_R[1][i_elec][i_spin]) 
                            W2 = tt.kron(II_fermibath, self.d_sys_L[0][i_elec][i_spin]) - tt.kron(II_fermibath_sign, self.d_sys_R[0][i_elec][i_spin]) 
                            XOP += -elec_lead_couplings[i_elec][i_spin][i_lead]**2/4.* ( W1*W2 + W2*W1)
                            XOP = XOP.round(self.eps)
    
                            OpMat = tt.kron(II_fermibath, self.d_sys_L[0][i_elec][i_spin] * self.d_sys_L[1][i_elec][i_spin])
                            Current[i_cstate] += -elec_lead_couplings[i_elec, i_spin, i_lead]**2/2.*tt.matvec(OpMat, Ix0)
                            OpMat = tt.kron(II_fermibath, self.d_sys_L[1][i_elec][i_spin] * self.d_sys_L[0][i_elec][i_spin])
                            Current[i_cstate] += elec_lead_couplings[i_elec, i_spin, i_lead]**2/2.*tt.matvec(OpMat, Ix0)
                        for i_pole in range(self.n_pole_fermi):
                            for i_sigma in range(2):
                                i_pos += 1
    
                                BB = op.OP1D_tt(self.mat_occupation, I_fermibath_list, i_pos, N_env)
                                Gop += -exponents_fermi[i_lead][i_spin][i_pole][i_sigma] * tt.kron(BB, self.II_sys)
    
                                Bp = op.OP1D_tt_fermi(self.mat_annihilation, I_fermibath_sign_list, I_fermibath_list, i_pos, N_env)
                                Sop += -1j * elec_lead_couplings[i_elec, i_spin, i_lead] * tt.kron(Bp, self.d_sys_L[1-i_sigma][i_elec][i_spin])
                                Bp = op.OP1D_tt_fermi(self.mat_annihilation, I_fermibath_list, I_fermibath_sign_list, i_pos, N_env)
                                Sop += 1j * elec_lead_couplings[i_elec, i_spin, i_lead] * tt.kron(Bp, self.d_sys_R[1-i_sigma][i_elec][i_spin])
                                Bm = op.OP1D_tt_fermi(self.mat_creation, I_fermibath_sign_list, I_fermibath_list, i_pos, N_env)
                                Rop += -1j * elec_lead_couplings[i_elec, i_spin, i_lead] * coefficients_fermi[i_lead][i_spin][i_pole] * tt.kron(Bm, self.d_sys_L[i_sigma][i_elec][i_spin])
                                Bm = op.OP1D_tt_fermi(self.mat_creation, I_fermibath_list, I_fermibath_sign_list, i_pos, N_env)
                                Rop += 1j * elec_lead_couplings[i_elec, i_spin, i_lead] * np.conjugate( coefficients_fermi[i_lead][i_spin][i_pole] ) * tt.kron(Bm, self.d_sys_R[i_sigma][i_elec][i_spin])
                                OpMat = tt.kron(II_fermibath, self.d_sys_L[i_sigma][i_elec][i_spin]) 
                                Ix=tt.mkron(v1[i_pos], self.I1_elec, self.I1_vib)
                                Current[i_cstate] += -1j* (-1)**i_sigma * elec_lead_couplings[i_elec, i_spin, i_lead] * tt.matvec(OpMat, Ix)
    
                                Gop = Gop.round(self.eps)
                                Sop = Sop.round(self.eps)
                                Rop = Rop.round(self.eps)
                        Current[i_cstate] = Current[i_cstate].round(self.eps)
        XOP += Gop
        XOP += Sop
        XOP += Rop
        XOP = XOP.round(self.eps)
       
        return XOP, Current

    def HQME_fermi_env_imagtime(self, coefficients_imag, exponents_imag, elec_lead_couplings, tau):

        XOP_imagtime = None
        XOP_imagtime += -tt.kron(self.II_fermibath_imagtime, self.Liouville_L)
        XOP_imagtime = XOP_imagtime.round(self.eps)
    
        i_pos = -1
        Sop = None
        Rop = None
        if self.HEAD_TAIL_COUP:
            for i_spin in range(self.n_spin):
                for i_lead in range(self.n_lead):
                    i_state = i_spin*self.n_lead+i_lead
                    for i_pole in range(self.n_pole_fermi):
                        for i_sigma in range(2):
                            for i_sincos in range(2):
                                i_pos += 1
                                if i_pole == 0:
                                    coeff_1 = -elec_lead_couplings[i_spin][i_lead] * coefficients_imag[i_lead][i_spin][i_pole] * ( (1-i_sincos)*np.sin(exponents_imag[i_lead][i_spin][i_pole]*(self.beta/2-tau)) - i_sincos*np.cos(exponents_imag[i_lead][i_spin][i_pole]*(self.beta/2-tau)) )
                                else:
                                    coeff_1 = -elec_lead_couplings[i_spin][i_lead] * coefficients_imag[i_lead][i_spin][i_pole] * ( (1-i_sincos)*np.cos(exponents_imag[i_lead][i_spin][i_pole]*tau) - i_sincos*np.sin(exponents_imag[i_lead][i_spin][i_pole]*tau) )
                                coeff_2 = elec_lead_couplings[i_spin][i_lead] * ( (1-i_sincos)*np.sin(exponents_imag[i_lead][i_spin][i_pole]*tau) + i_sincos*np.cos(exponents_imag[i_lead][i_spin][i_pole]*tau) )
    
                                Sop += coeff_1 * tt.kron(self.Bp_imagtime[i_pos], self.d_sys_L[1-i_sigma][i_state])
                                Rop += coeff_2 * tt.kron(self.Bm_imagtime[i_pos], self.d_sys_L[i_sigma][i_state])
                                Sop = Sop.round(self.eps)
                                Rop = Rop.round(self.eps)
        else:
            for i_elec in range(self.n_elec):
                for i_spin in range(self.n_spin):
                    for i_lead in range(self.n_lead):
                        for i_pole in range(self.n_pole_fermi):
                            for i_sigma in range(2):
                                for i_sincos in range(2):
                                    i_pos += 1
                                    if i_pole == 0:
                                        coeff_1 = -elec_lead_couplings[i_elec, i_spin, i_lead] * coefficients_imag[i_lead][i_spin][i_pole] * ( (1-i_sincos)*np.sin(exponents_imag[i_lead][i_spin][i_pole]*(self.beta/2-tau)) - i_sincos*np.cos(exponents_imag[i_lead][i_spin][i_pole]*(self.beta/2-tau)) )
                                    else:
                                        coeff_1 = -elec_lead_couplings[i_elec, i_spin, i_lead] * coefficients_imag[i_lead][i_spin][i_pole] * ( (1-i_sincos)*np.cos(exponents_imag[i_lead][i_spin][i_pole]*tau) - i_sincos*np.sin(exponents_imag[i_lead][i_spin][i_pole]*tau) )
                                    coeff_2 = elec_lead_couplings[i_elec, i_spin, i_lead] * ( (1-i_sincos)*np.sin(exponents_imag[i_lead][i_spin][i_pole]*tau) + i_sincos*np.cos(exponents_imag[i_lead][i_spin][i_pole]*tau) )
    
                                    Sop += coeff_1 * tt.kron(self.Bp_imagtime[i_pos], self.d_sys_L[1-i_sigma][i_elec][i_spin])
                                    Rop += coeff_2 * tt.kron(self.Bm_imagtime[i_pos], self.d_sys_L[i_sigma][i_elec][i_spin])
                                    Sop = Sop.round(self.eps)
                                    Rop = Rop.round(self.eps)
       
        XOP_imagtime += Sop
        XOP_imagtime += Rop
        XOP_imagtime = XOP_imagtime.round(self.eps)
    
        return XOP_imagtime

    

