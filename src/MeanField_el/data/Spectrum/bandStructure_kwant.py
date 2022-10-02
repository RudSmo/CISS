import kwant
import scipy.sparse.linalg as sla
import tinyarray

s0 = tinyarray.array([[1,0],[0,1]])


t0 = 1.0
tso = 1.0
ahel = 1.0
chel = 1.0
num_laps = 10
ions = 4
n_sites = int(num_laps*ions)

def UnitVector(a,c,j,num_laps,N):
        phij = 2*np.pi*(j-1)*num_laps/(N-1)
        cj = c*(j-1)/(N-1)
        return np.array([a*np.cos(phij),a*np.sin(phij),cj])

def VDot(a):
        sx = np.array([[0,1],[1,0]])
        sx = np.array([[0,-1j],[1j,0]])
        sz = np.array([[1,0],[0,-1]])
        return a[0]*sx + a[1] * sy + a[2] * sz

def SOC(site1,site2,tso,ahel,chel,num_laps,N):
	j = site1.pos[0]
	rj = UnitVector(ahel,chel,j,num_laps,N)
	rj1 = UnitVector(ahel,chel,j+1,num_laps,N)
	rj2 = UnitVector(ahel,chel,j+2,num_laps,N)
	rjrj1 = rj-rj1
	rjrj2 = rj-rj2
	dj1 = rjrj1/np.linalg.norm(rjrj1)
	dj2 = rjrj2/np.linalg.norm(rjrj2)
	vj = np.cross(dj1,dj2)
	Vj = 1j*tso*VDot(vj)
	return Vj

	


a = 1
lat = kwant.lattice.square()
syst = kwant.Builder()
syst[kwant.builder.HoppingKind((1,0), lat, lat)] = -t0*s0
syst[kwant.builder.HoppingKind((2,0), lat, lat)] = SOC
kwant.plot(syst.finalized())
#syst.finalized()
#ham_mat = syst.hamiltonian_submatrix(params=dict(tso=tso,ahel=ahel,chel=chel,num_laps=num_laps,N=n_sites),sparse=True)
#ev = sla.eigsh(ham_mat.tocsc(), k=15, sigma=0,return_eigenvectors=False)
#plt.plot(ev)
#plt.show()
