module MultiOpenIntegrator
include("global_constants.jl")
include("hamiltonian.jl")
include("iosave.jl")
include("electrons.jl")
include("rungekutta.jl")
include("classical_phonons.jl")
include("structures.jl")
using .GlobalConstants
using .Hamiltonian
using .IOsave
using .Electrons
using .RungeKutta
using .ClassicalPhonons
using .CustomStructures
using LinearAlgebra
using HDF5

export Propagate


function Propagate(
		  HamPar 	::Any,  
		  dt		::Number,
		  numstep	::Int,
		  tL		::Array,
		  tR		::Array,
		  beta		::Float64,
		  nPoles	::Int,
		  V         ::Float64,
		  path		::String,
		  rhoIn     ::String,
		  phIn	    ::String,
		  vSampl    ::String,
		  popIn 	= zeros(1),
		  pathIn    = "")


	if phIn == "TDSCF"
		X0 		= HamPar.X
		P0 		= HamPar.P
		weight  = prod(Wignerweights(X0,P0,beta,HamPar.ws,vSampl))
	elseif phIn == "Wigner"
		#A       = [SampleWigner(beta,HamPar.ws[i],vSampl) for i in 1:HamPar.N]
		#X0 		= [A[i][1] for i in 1:HamPar.N]
		#P0 		= [A[i][2] for i in 1:HamPar.N]

		X0,P0   = SampleWigner_Multi(beta,HamPar.ws,HamPar.N)
		weight  = prod(Wignerweights(X0,P0,beta,HamPar.ws,vSampl))
		
	elseif phIn == "fromFile"
		fid     = h5open(pathIn,"r")
		X0      = read(fid["q"])
		P0      = read(fid["p"])
		weight  = prod(Wignerweights(X0,P0,beta,HamPar.ws,vSampl))
	end

	Xin         = X0
	Pin         = P0
	X_l 		= [Xin]
	P_l 		= [Pin]
	HamPar.X    = Xin
	HamPar.P    = Pin

	hcen_i      = HamPar.Hamiltonian(HamPar)
	Gamma1,Gamma2 = GammaMatrices1d(HamPar.N,tL,tR)
	GammaT		= Gamma1+Gamma2
	tim   		= [0.0]
	kappa,xi	= PadeFermi(nPoles)
	chi1		= 0.0 .+ 1im .* xi ./beta    		# At t<t_0, voltage drop is zero.
	chi2		= 0.0 .+ 1im .* xi ./beta
	chiLP 	    = [chi1,chi2]
	if rhoIn == "pure"
		rho0    = Diagonal(float(popIn))
	elseif rhoIn == "equilibrium"
		rho0	= DensityFExOpen(HamPar.N,1,hcen_i,kappa,chiLP,GammaT,beta)*(1+0*1im)
	elseif rhoIn == "closedEquilibrium"
		rho0    = DensityFEx(HamPar.N,1,hcen_i,HamPar.N,0)
	elseif rhoIn == "fromFile"
		rho0    = read(fid["Rho"])
	end

	pop_l       = [Populations(rho0)]
    rho_i       = rho0
    
	
	if rhoIn == "fromFile"
		curp1_i = [read(fid["Curp1$(j)"]) for j in 1:nPoles]
		curp2_i = [read(fid["Curp2$(j)"]) for j in 1:nPoles]
		close(fid)
	else
		curp1_i = [CurrentInit_PFree(hcen_i,GammaT,chi1[i],kappa[i],beta,Gamma1) for i in 1:length(chi1)]
		curp2_i = [CurrentInit_PFree(hcen_i,GammaT,chi2[i],kappa[i],beta,Gamma2) for i in 1:length(chi2)]
	end

	curM1		= CurrentMatrix(rho0,Gamma1,curp1_i)
	curM2		= CurrentMatrix(rho0,Gamma2,curp2_i)
	cur1_l		= [Current(curM1)]
 	cur2_l		= [Current(curM2)]
	scur1_l 	= [SpinCurrent(curM1)]
	scur2_l 	= [SpinCurrent(curM2)]
	mom_el      = [MomentumElectrons(rho0)]
	loc_cur_l   = [LocalCurrent(rho0,hcen_i)]

	chi1	= V/2 .+ 1im .* xi ./beta 						# Overwrite previous chis by the new one WITH voltage drop
	chi2	= -V/2 .+ 1im .* xi ./beta
	
	for i in 1:numstep
		curp1_c, curp2_c, cur1_c, cur2_c, cur_c = RK_Currents(beta,hcen_i,Gamma1,Gamma2,GammaT,kappa,chi1,chi2,curp1_i,curp2_i,rho_i,dt)
		rho_c	= RK_VneumannL(rho_i,hcen_i,cur_c,dt)

		Xn,Pn   = SymplecticEulerStep(HamPar,rho_c,dt) 
		# BUFFER
		rho_i   = copy(rho_c)
		curp1_i = copy(curp1_c)
		curp2_i = copy(curp2_c)
		Xin     = Xn
		Pin     = Pn
		HamPar.X = Xn
		HamPar.P = Pn
		hcen_i  = HamPar.Hamiltonian(HamPar)	
		pop_l,tim,cur1_l,cur2_l,scur1_l,scur2_l,X_l,P_l,mom_el,loc_cur_l = SaveOutput(cur1_c,cur2_c,pop_l,rho_i,hcen_i,Xn,Pn,tim,i,dt,cur1_l,cur2_l,scur1_l,scur2_l,loc_cur_l,X_l,P_l,mom_el,weight,path)
		if i == numstep 
			fil = h5open("$(path)_final_Quantities","w")
			fil["Rho"] = rho_i
			for j in 1:nPoles
				fil["Curp1$(j)"] = curp1_i[j]
				fil["Curp2$(j)"] = curp2_i[j]
			end
			Xin = float.(Xin)
			Pin = float.(Pin)
			fil["X"] = Xin
			fil["P"] = Pin
			close(fil)	 
		end
	end
end

end
