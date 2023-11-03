
import numpy as np
from src import mashf90
from scipy.linalg import expm


class rpmd():
    def __init__(self,mass,beta,nb):
        self.beta=beta
        self.nb = nb
        self.betaN = float(beta)/nb
        self.mass = mass
        self.nf = len(mass)
    
    def step(self,p,q,dt):
        # TODO: free RP evolution
        p = p + self.force * dt/2
        q = q + p/self.mass * dt
        self.calcV(q)
        p = p + self.force * dt/2
        return p,q

    def anderson(self,p):
        psig = np.sqrt(self.mass/self.betaN)
        for k in range(self.nb):
            p[k] = np.random.normal(0,psig)

class MFrpmd(rpmd):
    def __init__(self,mass,beta,nb,ns):
        rpmd.__init__(self,mass,beta,nb)
        self.ns = ns

    def runtrj(self,dt,nt):
        """ Run equilibrium trajectory """
        nb = self.nb
        nf = self.nf
        ns = self.ns
        q = np.zeros((nb,nf))
        p = np.zeros((nb,nf))
        self.anderson(p)
        self.calcV(q)
        # Vt = np.zeros(nt+1)
        Et = np.zeros(nt+1)
        Pt = np.zeros((nt+1,ns))
        # qt = np.zeros((nt+1,nb,nf))
        # pt = np.zeros((nt+1,nb,nf))
        # Vt[0] = self.pot
        Et[0] = self.calcE(p)
        Pt[0] = self.calcP()
        # qt[0] = q.copy()
        # pt[0] = p.copy()
        for it in range(nt):
            p,q = self.step(p,q,dt)
            # TODO: Thermostat
            if it%int(np.sqrt(nt))==0:
                self.anderson(p)
            # print(q[0,:5])
            # qt[it+1] = q.copy()
            # pt[it+1] = p.copy()
            # Vt[it+1] = self.pot
            Et[it+1] = self.calcE(p)
            Pt[it+1] = self.calcP()
        return Et, Pt
    
    def calcE(self,p):
        return self.pot + 0.5*np.sum(np.sum(p**2,0)/self.mass)

    def calcV(self,q):
        nb = self.nb
        ns = self.ns
        nf = self.nf
        V = np.zeros((nb,ns,ns))
        dVdq = np.zeros((nb,nf,ns,ns))
        for k in range(nb):
            V[k],dVdq[k] = mashf90.pot_mat(q[k],ns)
        M = self.calcM(V)
        self.M = M
        self.pot = -1./self.betaN * np.log(self.calcZ(M))
        self.force = self.calcForce(M,dVdq)
        # """ Test derivative """
        # qtmp = q.copy()
        # dq = 1e-6
        # qtmp[0,0] = q[0,0] + dq
        # for k in range(nb):
        #     V[k],dVdq[k] = mashf90.pot_mat(qtmp[k],ns)
        # M2 = self.calcM(V)
        # pot2 = -1./self.betaN * np.log(self.calcZ(M2))
        # dpot = (pot2-self.pot)/dq
        # print(self.force[0,0],-dpot)

    def calcP(self):
        P = np.zeros(self.ns)
        for n in range(self.ns):
            nn = np.zeros((self.ns,self.ns))
            nn[n,n] = 1.
            for k in range(self.nb):
                nn = np.dot(nn,self.M[k])
            P[n] += np.trace(nn)
        P = P / self.calcZ(self.M)
        return P

    def calcM(self,V):
        return np.array([expm(-self.betaN*V[k]) for k in range(self.nb)])
        # M = np.zeros((self.nb,self.ns,self.ns))
        # for i in range(self.nb):
        #     M[i] = expm(-self.betaN*V[i])
        # return M

    def calcZ(self,M):
        rho = np.eye(self.ns)
        for i in range(self.nb):
            rho = np.dot(rho,M[i])
        return np.trace(rho)

    def calcForce(self,M,dVdq):
        """ Calculate tr[D_k H_k] where 
        D_k = dM_k/dq_k is the derivative of M_k=exp{-bN H(q_k)} and
        H_k = M_{k+1} ... M_n M_1 ... M_{k-1} is the 'hole' matrix.
        """
        D = -self.betaN*np.einsum('kfnm,kml->kfnl',dVdq,M)
        H = self.calcH(M)
        Z = self.calcZ(M)
        return np.einsum('kfnm,kmn->kf',D,H)/(self.betaN*Z)

    def calcH(self,M):
        """ Use Bell's algorithm to calculate 'hole' matrix from
        H_k = G_{k+1}F_{k-1}, which is built on the noration
        G_k = M_k ... M_n (right product) and
        F_k = M_1 ... M_k (left product) 
        """
        nb = self.nb
        ns = self.ns
        H = np.zeros((nb,ns,ns))
        if nb==1:
            H[0] = np.eye(ns)
            return H
        F = np.zeros((nb,ns,ns))
        G = np.zeros((nb,ns,ns))
        F[0] = M[0]
        G[-1] = M[-1]
        for k in range(1,nb):
            F[k] = np.dot(F[k-1],M[k])
            G[-1-k] = np.dot(M[-1-k],G[-k])
        H[0] = G[1]
        H[-1] = F[-2]
        for k in range(1,nb-1):
            H[k] = np.dot(G[k+1],F[k-1])
        return H
    




