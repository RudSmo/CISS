# Tutorial 2.2.2. Transport through a quantum wire
# ================================================
#
# Physics background
# ------------------
#  Conductance of a quantum wire; subbands
#
# Kwant features highlighted
# --------------------------
#  - Builder for setting up transport systems easily
#  - Making scattering region and leads
#  - Using the simple sparse solver for computing Landauer conductance

from matplotlib import pyplot
import kwant
import numpy as np
import tinyarray
# First, define the tight-binding system

sys = kwant.Builder()

# Here, we are only working with square lattices
a = 1
lat = kwant.lattice.square(a)

W = 1
t0 = 1.0
tso = 1.0
ahel = 1.0
chel = 1.0
num_laps = 10
ions = 4
L = int(num_laps*ions)
s0 = tinyarray.array([[1,0],[0,1]])


def UnitVector(a,c,j,num_laps,N):
	phij = 2*np.pi*(j-1)*num_laps/(N-1)
	cj = c*(j-1)/(N-1)
	return np.array([a*np.cos(phij),a*np.sin(phij),cj])

def VDot(a):
	sx = np.array([[0,1],[1,0]])
	sy = np.array([[0,-1j],[1j,0]])
	sz = np.array([[1,0],[0,-1]])
	return a[0]*sx + a[1] * sy + a[2] * sz

def SOC(j,tso,ahel,chel,num_laps,N):
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

# Define the scattering region

for i in range(L):
	for j in range(W):
		# On-site Hamiltonian
		sys[lat(i, j)] = 0.0

		# Hopping in y-direction
		if j > 0:
			sys[lat(i, j), lat(i, j - 1)] = 0.0

		# Hopping in x-direction
		if i > 0:
			sys[lat(i, j), lat(i - 1, j)] = -t0*s0

		if i < L-2:
			sys[lat(i, j), lat(i-2,j)] = SOC(i,tso,ahel,chel,num_laps,L)
	   
sys = sys.finalized()
kwant.plot(sys)
# Then, define and attach the leads:

# First the lead to the left
# (Note: TranslationalSymmetry takes a real-space vector)
#sym_left_lead = kwant.TranslationalSymmetry((-a, 0))
#left_lead = kwant.Builder(sym_left_lead)
#
#for j in range(W):
#    left_lead[lat(0, j)] = 4 * t
#    if j > 0:
#	left_lead[lat(0, j), lat(0, j - 1)] = -t
#    left_lead[lat(1, j), lat(0, j)] = -t
#
#sys.attach_lead(left_lead)
#
## Then the lead to the right
#sym_right_lead = kwant.TranslationalSymmetry((a, 0))
#right_lead = kwant.Builder(sym_right_lead)
#
#for j in range(W):
#    right_lead[lat(0, j)] = 4 * t
#    if j > 0:
#	right_lead[lat(0, j), lat(0, j - 1)] = -t
#    right_lead[lat(1, j), lat(0, j)] = -t

#sys.attach_lead(right_lead)

# Plot it, to make sure it's OK

# Finalize the system
#sys = sys.finalized()

# Now that we have the system, we can compute conductance
#energies = []
#data = []
#for ie in range(1000):
#    energy = ie * 0.01
#
#    # compute the scattering matrix at a given energy
#    smatrix = kwant.smatrix(sys, energy)
#
#    # compute the transmission probability from lead 0 to
#    # lead 1
#    energies.append(energy)
#    data.append(smatrix.transmission(1, 0))
#
## Use matplotlib to write output
## We should see conductance steps
#pyplot.figure()
#pyplot.plot(energies, data)
#pyplot.xlabel("energy [t]")
#pyplot.ylabel("conductance [e^2/h]")
#pyplot.show()
