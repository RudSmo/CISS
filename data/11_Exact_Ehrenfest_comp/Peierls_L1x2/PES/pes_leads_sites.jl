using LinearAlgebra
using PyPlot


function Hamiltonian(Nl,X,G0,lambda,e0,tL)
	# builds the entire hamiltonian for lead-sys-lead
	#     [ HL HLC HLR ]
	# H = | HCL HC HCR |
	#     [ HRL HRC HR ]
	# HLR = HRL = 0

	H  = zeros(Nl+2,Nl+2) 	# allocate Hamiltonian
	
	# Fill up HL with lead energies
	for i in 1:Nl-1
		H[i+1,i] = -tL 
		H[i,i+1] = -tL
	end
	
#	# Fill up HR with lead energies
#	for i in Nl+3:2*Nl+1
#		H[i+1,i] = -tL
#		H[i,i+1] = -tL
#	end

	# Fill up central Hamiltonian
	H[Nl+1,Nl+1] = e0 + X^2/2 # unoccup.
	H[Nl+2,Nl+2] = e0 + lambda*X + X^2/2 # occup.

	# Fill up HLC 
	H[Nl,Nl+1] = G0 * (1+exp(-X^2))
	H[Nl+1,Nl] = G0 * (1+exp(-X^2))	
	H[Nl,Nl+2] = G0 * (1+exp(-X^2))
	H[Nl+2,Nl] = G0 * (1+exp(-X^2))
	
#	H[Nl+3,Nl+2] = G0 * (1+exp(-X^2))
#	H[Nl+2,Nl+3] = G0 * (1+exp(-X^2))
#	H[Nl+3,Nl+1] = G0 * (1+exp(-X^2))
#	H[Nl+1,Nl+3] = G0 * (1+exp(-X^2))
	return H
end

Nl = 20
G0 = 0.0#3
e0 = 0.01875
lambda = sqrt(2) * 0.0075
Xrng = -1:0.001:1
Emin = 1.0

PES = zeros(Int(Nl+2),length(Xrng))
for (i,X) in enumerate(Xrng)
	H = Hamiltonian(Nl,X,G0,lambda,e0,Emin)
	PES[:,i] = eigvals(H)
end
for i in 1:Nl+2
	plt.plot(Xrng,PES[i,:])
end
plt.show()
