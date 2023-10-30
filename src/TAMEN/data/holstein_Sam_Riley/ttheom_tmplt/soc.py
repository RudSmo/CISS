import numpy as np

class spin_orbital_coupling():
    def __init__(self, a_helix, c_helix, number_laps, number_elec, chirality):
        self.a_helix = a_helix
        self.c_helix = c_helix
        self.number_laps = number_laps
        self.number_elec = number_elec
        self.chirality = chirality
    
    def VSI(self,n):
        def HelixVector(m):
            phi_m = 2*np.pi * (m-1) * self.number_laps / (self.number_elec - 1)
            return np.array([self.a_helix * np.cos( phi_m ), 
                             self.a_helix * np.sin( phi_m ),
                             self.c_helix * (m-1) / (self.number_elec - 1)])

        if self.chirality == "+":
            k1 = n+1
            k2 = n+2
        else:
            k1 = n-1
            k2 = n-2 
        r_n   = HelixVector(n)
        r_n1  = HelixVector(k1)
        r_n2  = HelixVector(k2)

        d_n1u = r_n - r_n1
        d_n1  = d_n1u / np.sqrt( np.dot( d_n1u,d_n1u ) )
        d_n2u = r_n - r_n2
        d_n2  = d_n2u / np.sqrt( np.dot( d_n2u,d_n2u ) )

        v_n   = np.cross( d_n1, d_n2 )
        sx    = np.array([[0, 1.],[1., 0]])
        sy    = np.array([[0,-1j],[1j, 0]])
        sz    = np.array([[1, 0.],[0.,-1]])
        return v_n[0] * sx + v_n[1] * sy + v_n[2] * sz 
     
if __name__ == "__main__":
    chirality = "+"
    r_helix = 1.0 
    l_helix = 1.0
    M_laps = 2
    Dim_elecDOF = M_laps*2
    so_coup = spin_orbital_coupling(r_helix, l_helix, M_laps, Dim_elecDOF, chirality)
    print(so_coup.VSI(1))
