module MultiOpenIntegrator
include("global_constants.jl")
include("hamiltonian.jl")
include("iosave.jl")
include("electrons.jl")
include("rungekutta.jl")
include("classical_phonons.jl")
include("structures.jl")
include("observables.jl")
include("pade.jl")
using .GlobalConstants
using .Hamiltonian
using .IOsave
using .Electrons
using .RungeKutta
using .ClassicalPhonons
using .CustomStructures
using .Observables
using .Pade
using LinearAlgebra
using HDF5

export Propagate


function Propagate(
		  HamPar 	::Any,  
		  dt		::Number,
		  numstep	::Int,
		  tL		::Array,
		  tR		::Array,
		  BetaL		::Array,
		  DeltaL	::Array,
		  nPoles	::Int,
		  Vs         	::Tuple,
		  Vtype		::String,
		  Vdrop 	::String,
		  path		::String,
		  rhoIn     	::String,
		  phIn	    	::String,
		  vSampl    	::String,
		  popIn 	= zeros(1),
		  pathIn    	= "")

	# Quantities:
	#
	# Reduced single-particle density matrix
	# rho = Array(nSys,nSys)
	#
	# Current matrices
	# Pi  = Array(nLeads,nSys,nSys)  
	#
	# Auxiliary mode vectors
	# Psi = Array(nLeads,nPoles,nChannels,nSys) 
	#
	# Channel-resolved coupling / Eigenvalues of coupling matrices
	# Gamma = Array(nLeads,nChannels)
	#
	# Eigenvectors of coupling matrices
	# Xi  = Array(nLeads,nChannels,nSys) 
	#
	# Auxiliary two-mode coefficients
	# Omega = Array(nLeads,nPoles,nChannels,nLeads,nPoles,nChannels)
	#
	# Pade poles
	# ChiPade = Array(nLeads,nPoles,nSign) with nSign = 1 -> "µ+iz_p/beta_a" and nSign = 2 -> "µ-iz_p/beta_a"
	#
	# Pade coefficients
	# EtaPade = Array(nPoles)
	# 
	# Inverse temperature of leads
	# BetaL = Array(nLeads)
	#
	# Electronic energy shift of leads
	# DeltaL = Array(nLeads,numstep)

	nSys = 2*HamPar.N

	# This code is assuming two leads, one left and one right
	nLeads = 2

	# The central system is coupled at the edges to the leads, and since there is spin-up and spin-down, 
	# there are two channels per lead
#	nChannels = nSys
	nChannels = 2

	nAll = [nSys,nLeads,nPoles,nChannels] 

	if phIn == "TDSCF"
		X0 		= HamPar.X
		P0 		= HamPar.P
	elseif phIn == "Wigner"
		#A      	= [SampleWigner(beta,HamPar.ws[i],vSampl) for i in 1:HamPar.N]
		#X0 		= [A[i][1] for i in 1:HamPar.N]
		#P0 		= [A[i][2] for i in 1:HamPar.N]

		X0,P0   	= SampleWigner_Multi(beta,HamPar.ws,HamPar.N)
		
	elseif phIn == "fromFile"
		fid     = h5open(pathIn,"r")
		X0      = read(fid["q"])
		P0      = read(fid["p"])
	end

	if Vtype == "const"
		V = Vs[1]*ones(numstep)
	elseif Vtype == "lorentz"
		TI = [i*dt for i in 1:numstep]
		t0,tau = Vs
		V = 2*tau ./ ((TI .- t0) .^ 2 .+ tau^2)	
	end
	µPot 	= zeros(nLeads,numstep)
	if Vdrop == "sym"
		µPot[1,:] = V ./ 2
		µPot[2,:] = - V ./ 2
	elseif Vdrop == "left"
		µPot[1,:] = V
		µPot[2,:] = zeros(numstep)
	end

	Xin         	= X0
	Pin         	= P0
	X_l 	    	= [Xin]
	P_l 	    	= [Pin]
	HamPar.X    	= Xin
	HamPar.P    	= Pin

	hcen_i      	= HamPar.Hamiltonian(HamPar)
	Gamma1,Gamma2 	= GammaMatrices1d(HamPar.N,tL,tR)
	GammaT		= Gamma1+Gamma2
	Gamma		= zeros(nLeads,nChannels)
#	Gamma[1,:]	= eigvals(Gamma1)
#	Gamma[2,:]	= eigvals(Gamma2)
	Gamma[1,1]	= tL[1]
	Gamma[1,2]	= tL[2]
	Gamma[2,1]	= tR[1]
	Gamma[2,2]	= tR[2]

	
	Axi		= zeros(nLeads,nChannels,nSys)
#	Axi[1,:,:]	= eigvecs(Gamma1)
#	Axi[2,:,:]	= eigvecs(Gamma2)
	Axi[1,1,1]      = 1.0	
	Axi[1,2,2] 	= 1.0
	Axi[2,1,end-1]  = 1.0
	Axi[2,2,end] 	= 1.0

	tim   		= [0.0]
	EtaPade,xi	= PadeFermi(nPoles)
	ChiPade		= zeros(nLeads,nPoles,2)*(1+0.0*1im)
	for a in 1:nLeads
		for p in 1:nPoles
			for s in 1:2
				ChiPade[a,p,s] = 0.0 + (-1)^(s-1) * 1im * xi[p]/BetaL[a] 
			end
		end
	end
			

	if rhoIn == "pure"
		rho0    = Diagonal(float(popIn))
	elseif rhoIn == "equilibrium"
		rho0	= Rho0_Equilibrium_Open(HamPar.N,hcen_i,nPoles,EtaPade,ChiPade,GammaT,BetaL)*(1+0*1im)
	elseif rhoIn == "closedEquilibrium"
		rho0    = Rho0_Equilibrium_Closed(HamPar.N,hcen_i,HamPar.N,0)
	elseif rhoIn == "fromFile"
		rho0    = read(fid["Rho"])
	end

	pop_l       = [Populations(rho0)]
	pop_es_l    = [PopulationsES(rho0,hcen_i)]
	rho_i       = rho0
   
	omega_i = zeros(nLeads,nPoles,nChannels,nLeads,nPoles,nChannels) * (1.0 + 0*1im) 
	psi_i 	= zeros(nLeads,nPoles,nChannels,nSys) * (1.0 + 0*1im)
	cur_i 	= CurrentMatrix(nAll,psi_i,Axi)	
		
	cur1_l		= [Current(cur_i[1,:,:])]
 	cur2_l		= [Current(cur_i[2,:,:])]
	scur1_l 	= [SpinCurrent(cur_i[1,:,:])]
	scur2_l 	= [SpinCurrent(cur_i[2,:,:])]
	mom_el      	= [MomentumElectrons(rho_i)]
	loc_cur_l   	= [LocalCurrent(rho_i,hcen_i)]
	
	for i in 1:numstep
		ChiPade = zeros(nLeads,nPoles,2) * (1+0.0*1im)
	        for a in 1:nLeads
        	        for p in 1:nPoles
                	        for s in 1:2
                        	        ChiPade[a,p,s] = µPot[a,i] + (-1)^(s-1) * 1im * xi[p]/BetaL[a]
	                        end
        	        end
	        end	
		omega_c, psi_c, rho_c = RK_Auxiliary(nAll,hcen_i,ChiPade,EtaPade,BetaL,DeltaL,Gamma,Axi,omega_i,psi_i,rho_i,dt)
		Xn,Pn    = SymplecticEulerStep(HamPar,rho_c,dt) 
		# BUFFER
	
		rho_i    = copy(rho_c)
		psi_i	 = copy(psi_c)
		omega_i	 = copy(omega_c)

		Xin      = Xn
		Pin      = Pn
		HamPar.X = Xn
		HamPar.P = Pn

		hcen_i   = HamPar.Hamiltonian(HamPar)	

		pop_l,pop_es_l,tim,cur1_l,cur2_l,scur1_l,scur2_l,X_l,P_l,mom_el,loc_cur_l = SaveOutput(nAll,psi_c,Axi,pop_l,pop_es_l,rho_i,hcen_i,Xn,Pn,tim,i,dt,cur1_l,cur2_l,scur1_l,scur2_l,loc_cur_l,X_l,P_l,mom_el,path)
		if i == numstep 
			fil = h5open("$(path)_final_Quantities","w")
			fil["Rho"] = rho_i
			fil["Omega"] = omega_i
			fil["Psi"] = psi_i
			Xin = float.(Xin)
			Pin = float.(Pin)
			fil["X"] = Xin
			fil["P"] = Pin
			close(fil)	 
		end
	end
end

end
