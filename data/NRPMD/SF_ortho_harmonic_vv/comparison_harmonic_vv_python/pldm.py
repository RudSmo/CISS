from __future__ import division

"""Built-in python imports"""
import numpy as np
import scipy.integrate

"""System import is achieved either by
	- pldm.py being run from directory containing system.py
	- the folder containing system.py being on $PYTHPATH"""
import system

"""Local imports"""
from analysis import Observables
from arguments_pldm import arguments
import model_pldm

"""Lib-python imports"""
from polylib.PES import rp

"""f2py imports"""
from src import nrpmdf90

class Initialize(object):
	"""Initialize the necessary arrays for the traject object"""
	def __init__(self, args, model, sgn, polymer):
		"""Everything important happens here"""
		self.output = np.zeros((args.TS+1,4*model.S+2*model.F))
		self.output1 = np.zeros((args.TS+1,model.S,model.S),dtype=complex)
		self.output2 = np.zeros((args.TS+1,model.S,model.S),dtype=complex)
		self.sgn = sgn
		self.polymer = polymer

class Trajectory(object):
	"""Wrapper for all the data returned from FORTRAN trajectories"""
	def __init__(self, coord, coord1, coord2, args, model, j):
		"""Everything important happens here"""
		if (args.integrator == 'python'):
			"""Error tolerances for differential equation solver"""
			abserr = 1e-10
			relerr = 1e-10
			"""Place the initial Statonovich Weyl Kernels into the initial array"""
			coord_tot = np.zeros(4*model.S+2*model.F+4*(model.S**2))
			coord_tot[0:4*model.S+2*model.F] = coord
			coord_tot[4*model.S+2*model.F:4*model.S+2*model.F+(model.S**2)] = np.real(np.reshape(coord1,(1,model.S**2)))
			coord_tot[4*model.S+2*model.F+(model.S**2):4*model.S+2*model.F+2*(model.S**2)] = np.imag(np.reshape(coord1,(1,model.S**2)))
			coord_tot[4*model.S+2*model.F+2*(model.S**2):4*model.S+2*model.F+3*(model.S**2)] = np.real(np.reshape(coord2,(1,model.S**2)))
			coord_tot[4*model.S+2*model.F+3*(model.S**2):4*model.S+2*model.F+4*(model.S**2)] = np.imag(np.reshape(coord2,(1,model.S**2)))
			"""Perform the dynamics"""
			out = scipy.integrate.odeint(eom,coord_tot,t[j*args.TS2:(j+1)*args.TS2+1],args=(model.S,model.F,model.mass,model.zpe,),atol=abserr,rtol=relerr)
			"""Generate the time-evolved Stratonovich Weyl kernels"""
			self.output = out[:,0:4*model.S+2*model.F]
			self.output1 = 1j*np.reshape(out[:,4*model.S+2*model.F+(model.S**2):4*model.S+2*model.F+2*(model.S**2)],(args.TS2+1,model.S,model.S))
			self.output1 += np.reshape(out[:,4*model.S+2*model.F:4*model.S+2*model.F+(model.S**2)],(args.TS2+1,model.S,model.S)) 
			self.output2 = 1j*np.reshape(out[:,4*model.S+2*model.F+3*(model.S**2):4*model.S+2*model.F+4*(model.S**2)],(args.TS2+1,model.S,model.S))
			self.output2 += np.reshape(out[:,4*model.S+2*model.F+2*(model.S**2):4*model.S+2*model.F+3*(model.S**2)],(args.TS2+1,model.S,model.S)) 
		elif (args.integrator == 'harmonic'):
			self.output, self.output1, self.output2 = harmonic(coord,coord1,coord2,t[j*args.TS2:(j+1)*args.TS2+1],model.S,model.F,model.mass,model.omega,model.Vconst,model.Vlin,model.zpe)
		elif (args.integrator == 'harmonic_VV'):
			self.output, self.output1, self.output2 = harmonic_VV(coord,coord1,coord2,t[j*args.TS2:(j+1)*args.TS2+1],model.S,model.F,model.mass,model.omega,model.Vconst,model.Vlin,model.zpe)

def eom(x, t, S, F, mass, gamma):
	gradient = np.zeros(2*F+4*S+4*(S**2))

	"""Evolve the electronic mapping variables"""
	ham = nrpmdf90.pot(x[4*S:4*S+F], F, S)[1]
	ham = ham - (np.sum(np.diag(ham))/float(S))*np.identity(S)
	gradient[0:S] = np.dot(ham,x[S:2*S])
	gradient[S:2*S] = -np.dot(ham,x[0:S])
	gradient[2*S:3*S] = np.dot(ham,x[3*S:4*S])
	gradient[3*S:4*S] = -np.dot(ham,x[2*S:3*S])

	"""Evolve the Stratonovich-Weyl kernels"""
	coord1_re = np.reshape(x[4*S+2*F:4*S+2*F+(S**2)],(S,S)) 
	coord1_im = np.reshape(x[4*S+2*F+(S**2):4*S+2*F+2*(S**2)],(S,S))
	coord2_re = np.reshape(x[4*S+2*F+2*(S**2):4*S+2*F+3*(S**2)],(S,S))    
	coord2_im = np.reshape(x[4*S+2*F+3*(S**2):4*S+2*F+4*(S**2)],(S,S))
	output1_re = np.dot(ham,coord1_im)
	output1_im = -np.dot(ham,coord1_re)
	output2_re = np.dot(ham,coord2_im)
	output2_im = -np.dot(ham,coord2_re)
	gradient[4*S+2*F:4*S+2*F+(S**2)] = np.reshape(output1_re,(1,S**2))
	gradient[4*S+2*F+(S**2):4*S+2*F+2*(S**2)] = np.reshape(output1_im,(1,S**2))
	gradient[4*S+2*F+2*(S**2):4*S+2*F+3*(S**2)] = np.reshape(output2_re,(1,S**2))
	gradient[4*S+2*F+3*(S**2):4*S+2*F+4*(S**2)] = np.reshape(output2_im,(1,S**2))

	"""Evolve the nuclear degrees of freedom"""
	force0 = nrpmdf90.pot(x[4*S:4*S+F], F, S)[2]
	force1 = nrpmdf90.pot(x[4*S:4*S+F], F, S)[3]
	force2 = (np.einsum('ijj',force1)/float(S))*(1 - 0.25*np.sum(x[0:4*S]**2))
	force2 += 0.25*np.real(np.einsum('j,ij',x[0:S]+1j*x[S:2*S],np.einsum('j,ijk',x[0:S]-1j*x[S:2*S],force1))) 
	force2 += 0.25*np.real(np.einsum('j,ij',x[2*S:3*S]+1j*x[3*S:4*S],np.einsum('j,ijk',x[2*S:3*S]-1j*x[3*S:4*S],force1)))
	gradient[4*S:4*S+F] = x[4*S+F:4*S+2*F]/mass
	gradient[4*S+F:4*S+2*F] = -force0 - force2

	return gradient

def harmonic(x, x1, x2, t, S, F, mass, omega, Vconst, Vlin, gamma):

	"""Initialise arrays for time evolution"""
	tstep = len(t)
	dt = t[1] - t[0]
	output = np.zeros((tstep,4*S+2*F))
	output1 = np.zeros((tstep,S,S),dtype=complex)
	output2 = np.zeros((tstep,S,S),dtype=complex)
	output[0,:] = x
	output1[0,:,:] = x1
	output2[0,:,:] = x2

	coeff = np.zeros((S,2),dtype=complex)

	for i in range(0,tstep-1):
		"""Evolve the nuclear degrees of freedom exactly"""
		temp = np.copy(output[i,4*S:4*S+2*F])
		force = (np.einsum('ijj',Vlin)/float(S))*(1 - 0.25*np.sum(output[i,0:4*S]**2))
		force += 0.25*np.real(np.einsum('j,ij',output[i,0:S]+1j*output[i,S:2*S],np.einsum('j,ijk',output[i,0:S]-1j*output[i,S:2*S],Vlin))) 
		force += 0.25*np.real(np.einsum('j,ij',output[i,2*S:3*S]+1j*output[i,3*S:4*S],np.einsum('j,ijk',output[i,2*S:3*S]-1j*output[i,3*S:4*S],Vlin)))
		"""Evolve nuclear positions"""
		output[i+1,4*S:4*S+F] = (temp[0:F]+(force/(mass*(omega**2))))*np.cos(0.5*omega*dt) + (temp[F:2*F]/(mass*omega))*np.sin(0.5*omega*dt) - (force/(mass*(omega**2))) 
		"""Evolve nuclear momenta"""
		output[i+1,4*S+F:4*S+2*F] = -(temp[0:F]+(force/(mass*(omega**2))))*mass*omega*np.sin(0.5*omega*dt) + temp[F:2*F]*np.cos(0.5*omega*dt) 

		"""Generate the electronic Hamiltonian"""
		ham = Vconst + np.einsum('i,ijk',output[i+1,4*S:4*S+F],Vlin)
		ham = ham - (np.sum(np.diag(ham))/float(S))*np.identity(S)
		
		"""Evolve the electronic degrees of freedom exactly""" 
		prop = scipy.linalg.expm(-1j*ham*dt)
		coeff[:,0] = np.dot(prop,output[i,0:S]+1j*output[i,S:2*S])
		coeff[:,1] = np.dot(prop,output[i,2*S:3*S]+1j*output[i,3*S:4*S])
		output1[i+1,:,:] = np.matmul(prop,output1[i,:,:])
		output2[i+1,:,:] = np.matmul(prop,output2[i,:,:])
		output[i+1,0:S] = np.real(coeff[:,0])
		output[i+1,S:2*S] = np.imag(coeff[:,0])
		output[i+1,2*S:3*S] = np.real(coeff[:,1]) # Mapping X
		output[i+1,3*S:4*S] = np.imag(coeff[:,1]) # Mapping P

		"""Evolve the nuclear degrees of freedom exactly"""
		temp = np.copy(output[i+1,4*S:4*S+2*F])
		force = (np.einsum('ijj',Vlin)/float(S))*(1 - 0.25*np.sum(output[i+1,0:4*S]**2))
		force += 0.25*np.real(np.einsum('j,ij',output[i+1,0:S]+1j*output[i+1,S:2*S],np.einsum('j,ijk',output[i+1,0:S]-1j*output[i+1,S:2*S],Vlin))) 
		force += 0.25*np.real(np.einsum('j,ij',output[i+1,2*S:3*S]+1j*output[i+1,3*S:4*S],np.einsum('j,ijk',output[i+1,2*S:3*S]-1j*output[i+1,3*S:4*S],Vlin)))
		"""Evolve nuclear positions"""
		output[i+1,4*S:4*S+F] = (temp[0:F]+(force/(mass*(omega**2))))*np.cos(0.5*omega*dt) + (temp[F:2*F]/(mass*omega))*np.sin(0.5*omega*dt) - (force/(mass*(omega**2))) 
		"""Evolve nuclear momenta"""
		output[i+1,4*S+F:4*S+2*F] = -(temp[0:F]+(force/(mass*(omega**2))))*mass*omega*np.sin(0.5*omega*dt) + temp[F:2*F]*np.cos(0.5*omega*dt)  

	return output, output1, output2

def harmonic_VV(x, x1, x2, t, S, F, mass, omega, Vconst, Vlin, gamma):

	"""Initialise arrays for time evolution"""
	tstep = len(t)
	dt = t[1] - t[0]
	output = np.zeros((tstep,4*S+2*F))
	output1 = np.zeros((tstep,S,S),dtype=complex)
	output2 = np.zeros((tstep,S,S),dtype=complex)
	output[0,:] = x
	output1[0,:,:] = x1
	output2[0,:,:] = x2

	coeff = np.zeros((S,2),dtype=complex)

	for i in range(0,tstep-1):
		"""Evolve the nuclear degrees, first half step"""
		temp = np.copy(output[i,4*S:4*S+2*F])
		force = (np.einsum('ijj',Vlin)/float(S))*(1 - 0.25*np.sum(output[i,0:4*S]**2))
		force += 0.25*np.real(np.einsum('j,ij',output[i,0:S]+1j*output[i,S:2*S],np.einsum('j,ijk',output[i,0:S]-1j*output[i,S:2*S],Vlin)))
		force += 0.25*np.real(np.einsum('j,ij',output[i,2*S:3*S]+1j*output[i,3*S:4*S],np.einsum('j,ijk',output[i,2*S:3*S]-1j*output[i,3*S:4*S],Vlin)))
		"""Evolve nuclear momenta"""
		output[i+1,4*S+F:4*S+2*F] = temp[F:2*F] + 0.5 * force/mass * dt
		"""Evolve nuclear positions"""
		output[i+1,4*S:4*S+F] = temp[0:F] + output[i+1,4*S+F:4*S+2*F] * dt 

		"""Generate the electronic Hamiltonian"""
		ham = Vconst + np.einsum('i,ijk',output[i+1,4*S:4*S+F],Vlin)
		ham = ham - (np.sum(np.diag(ham))/float(S))*np.identity(S)

		"""Evolve the electronic degrees of freedom exactly"""
		prop = scipy.linalg.expm(-1j*ham*dt)
		coeff[:,0] = np.dot(prop,output[i,0:S]+1j*output[i,S:2*S])
		coeff[:,1] = np.dot(prop,output[i,2*S:3*S]+1j*output[i,3*S:4*S])
		output1[i+1,:,:] = np.matmul(prop,output1[i,:,:])
		output2[i+1,:,:] = np.matmul(prop,output2[i,:,:])
		output[i+1,0:S] = np.real(coeff[:,0])
		output[i+1,S:2*S] = np.imag(coeff[:,0])
		output[i+1,2*S:3*S] = np.real(coeff[:,1]) # Mapping X
		output[i+1,3*S:4*S] = np.imag(coeff[:,1]) # Mapping P

		"""Evolve the nuclear degrees of freedom exactly"""
		temp = np.copy(output[i+1,4*S:4*S+2*F])
		force = (np.einsum('ijj',Vlin)/float(S))*(1 - 0.25*np.sum(output[i+1,0:4*S]**2))
		force += 0.25*np.real(np.einsum('j,ij',output[i+1,0:S]+1j*output[i+1,S:2*S],np.einsum('j,ijk',output[i+1,0:S]-1j*output[i+1,S:2*S],Vlin)))
		force += 0.25*np.real(np.einsum('j,ij',output[i+1,2*S:3*S]+1j*output[i+1,3*S:4*S],np.einsum('j,ijk',output[i+1,2*S:3*S]-1j*output[i+1,3*S:4*S],Vlin)))
		"""Evolve nuclear momenta"""
		output[i+1,4*S+F:4*S+2*F] = temp[F:2*F] + 0.5 * force/mass * dt 

	return output, output1, output2


"""Read command line arguments"""
args = arguments()

"""Calculate the number of time steps for the total evolution"""
args.TS = (args.njump+1)*args.TS2

"""Load the system from system.py"""
model = system.Model(args)
model1 = model_pldm.Model(args, model)

"""State the norm associated with each Stratonovich-Weyl kernel"""
norms = {"standard-pldm": 1., "fullsphere" : model.S, "focused" : model.S}

"""Initialisation of ring polymer"""
print("\nRing Polymer init")
RP = rp.cart(args.N, None, model.mass, beta=args.beta)

"""Initialisation of system"""
print("System init")
model.initialize(args, RP.omega)

"""Set up data analysis"""
obs = Observables(args, model1)

"""Initial conditions for the trajectory"""
coord = np.zeros(4*model.S + 2*model.F)
t = np.linspace(0, args.TS*args.dt, args.TS+1)

"""numerical solution of the equation of motion differential equations"""
for i in range(0,args.traj): 

	"""Generate the first orthogonal coherent state basis"""
	coord[0:model.S], coord[model.S:2*model.S], sgn = model1.sample_electronic(args)
	"""Generate the Stratonovich-Weyl kernel"""
	coord1 = 0.5*(np.outer(coord[0:model.S]+1j*coord[model.S:2*model.S],coord[0:model.S]-1j*coord[model.S:2*model.S]) - model.zpe*np.identity(model.S))
	"""Generate the second orthogonal coherent state basis"""
	coord[2*model.S:3*model.S], coord[3*model.S:4*model.S], sgn = model1.sample_electronic(args)
	"""Generate the Stratonovich-Weyl kernel"""
	coord2 = 0.5*(np.outer(coord[2*model.S:3*model.S]+1j*coord[3*model.S:4*model.S],coord[2*model.S:3*model.S]-1j*coord[3*model.S:4*model.S]) - model.zpe*np.identity(model.S))
	"""Generate the initial nuclear coordinates"""
	coord[4*model.S:4*model.S+model.F], coord[4*model.S+model.F:4*model.S+2*model.F], polymer = model.sample_nuclear(args)

	"""Initialize the output arrays in the traject object"""
	traject = Initialize(args, model, sgn, polymer)

	"""Perform the classical trajectory dynamics for the first jump interval"""
	traject1=Trajectory(coord, coord1, coord2, args, model, 0)

	"""Calculate the output associated with the first jump interval"""
	traject.output[0:args.TS2+1,:] = traject1.output
	traject.output1[0:args.TS2+1,:,:] = norms[args.electronic_sampling]*traject1.output1
	traject.output2[0:args.TS2+1,:,:] = norms[args.electronic_sampling]*traject1.output2

	"""Perform jumps"""
	for j in range(1,args.njump+1):
		"""Generate the first orthogonal coherent state basis"""
		coord[0:model.S], coord[model.S:2*model.S], sgn = model1.sample_electronic(args)
		"""Generate the Startonovich Weyl kernel"""
		coord1 = 0.5*(np.outer(coord[0:model.S]+1j*coord[model.S:2*model.S],coord[0:model.S]-1j*coord[model.S:2*model.S]) - model.zpe*np.identity(model.S))
		"""Generate the second orthogonal coherent state basis"""
		coord[2*model.S:3*model.S], coord[3*model.S:4*model.S], sgn = model1.sample_electronic(args)
		"""Generate the Startonovich Weyl kernel"""
		coord2 = 0.5*(np.outer(coord[2*model.S:3*model.S]+1j*coord[3*model.S:4*model.S],coord[2*model.S:3*model.S]-1j*coord[3*model.S:4*model.S]) - model.zpe*np.identity(model.S))
		"""Set the initial nuclear coordinates to be those at the end of the last jump interval"""
		coord[4*model.S:4*model.S+2*model.F] = traject1.output[args.TS2,4*model.S:4*model.S+2*model.F]

		"""Perform the classical trajectory dynamics for the (j+1)th jump interval"""
		traject1=Trajectory(coord, coord1, coord2, args, model, j)

		"""Calculate the output associated with the (j+1)th jump interval"""
		traject.output[j*args.TS2:(j+1)*args.TS2+1,:] = traject1.output
		traject.output1[j*args.TS2:(j+1)*args.TS2+1,:,:] = norms[args.electronic_sampling]*np.dot(traject1.output1,traject.output1[j*args.TS2,:,:])
		traject.output2[j*args.TS2:(j+1)*args.TS2+1,:,:] = norms[args.electronic_sampling]*np.dot(traject1.output2,traject.output2[j*args.TS2,:,:])

	"""Calculate observables"""
	obs.add_traj(traject)

	"""Output to file"""
	if (i + 1) % args.iofreq == 0:
		print("Completed time steps: ", i+1)
		obs.restart()

"""Calculate the observable averages"""
obs.average()

"""Output to file"""
obs.output()
