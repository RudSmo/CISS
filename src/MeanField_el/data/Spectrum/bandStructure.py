from pythtb import *
import matplotlib.pyplot as plt
import numpy as np

lat = [[1.0]]
orb = [[0.0],[0.0]]
model = tb_model(1,1,lat,orb)
model.set_hop(-1.0,0,0,[1])
model.set_hop(-1.0,1,1,[1])
model.display()
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

def SOC(tso,ahel,chel,site,num_laps,N):
	rj = UnitVector(ahel,chel,site,num_laps,N)
	rj1 = UnitVector(ahel,chel,site+1,num_laps,N)
	rj2 = UnitVector(ahel,chel,site+2,num_laps,N)
	rjrj1 = rj-rj1
	rjrj2 = rj-rj2

	dj1 = rjrj1/np.linalg.norm(rjrj1)
	dj2 = rjrj2/np.linalg.norm(rjrj2)

	vj = np.cross(dj1,dj2)
	Vj = 1j*tso*VDot(vj)
	return Vj

sites = range(1,nsites+1)
for site in sites:
	tsoMat = SOC(tso,ahel,chel,site)
	model.set_hop(tsoMat[0,0],0,0,)
