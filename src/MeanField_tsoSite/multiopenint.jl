module MultiOpenIntegrator
include("global_constants.jl")
include("hamiltonian.jl")
include("iosave.jl")
include("electrons.jl")
include("rungekutta.jl")
include("classical_phonons.jl")
using .GlobalConstants
using .Hamiltonian
using .IOsave
using .Electrons
using .RungeKutta
using .ClassicalPhonons
using LinearAlgebra
using HDF5

export Propagate,PropagateSingle


function Propagate(N::Int,
		  E         ::Array,
		  t         ::Array,
		  tso       ::Array,
		  a         ::Float64,
		  c         ::Float64,
		  num_laps  ::Int64,
		  tph 		::Array,
		  tsoph 	::Array,  
		  ws 		::Array,
		  dt		::Number,
		  numstep	::Int,
		  tL		::Array,
		  tR		::Array,
		  beta		::Float64,
		  nPoles	::Int,
		  V         ::Array,
		  path		::String,
		  resampling::Bool,
		  rhoIn     ::String,
		  phIn	    ::String,
		  vSampl    ::String,
		  pathIn    = "",
		  popIn 	= zeros(1),
		  qnIn  	= zeros(1),
		  pnIn  	= zeros(1),
		  rsmpl_stp	= 0)


	if phIn == "TDSCF"
		q0 		= qnIn
		p0 		= pnIn
	elseif phIn == "Wigner"
		A       = [SampleWignerA(beta,ws[i],vSampl) for i in 1:N]
		q0 	    = [A[i][1] for i in 1:N]
		p0 		= [A[i][2] for i in 1:N]
	elseif phIn == "fromFile"
		fid     = h5open(pathIn,"r")
		q0      = read(fid["q"])
		p0      = read(fid["p"])
	end

	qin         = q0
	pin         = p0
	q_l 		= [q0]
	p_l 		= [p0]

	hcen_i      = FranssonHamiltonian(N,E,t,tso,a,c,num_laps,tph,tsoph,ws,q0)
	Gamma1,Gamma2 = GammaMatrices1d(N,tL,tR)
	GammaT		= Gamma1+Gamma2
	tim   		= [0.0]
	kappa,xi	= PadeFermi(nPoles)
	chi1		= V[1]/2 .+ 1im .* xi ./beta    		# At t<t_0, voltage drop is zero.
	chi2		= -V[1]/2 .+ 1im .* xi ./beta
	chiLP 	    = [chi1,chi2]
	if rhoIn == "pure"
		rho0    = Diagonal(float(popIn))
	elseif rhoIn == "equilibrium"
		rho0	= DensityFExOpen(N,1,hcen_i,kappa,chiLP,GammaT,beta)*(1+0*1im)
	elseif rhoIn == "closedEquilibrium"
		rho0    = DensityFEx(N,1,hcen_i,N,0)
	elseif rhoIn == "fromFile"
		rho0    = read(fid["Rho"])
	end

	pop_l       = [Populations(rho0)]
    rho_i       = rho0
    
	hcen0 		= hcen_i

	if rhoIn == "fromFile"
		curp1_i = [read(fid["Curp1$(j)"]) for j in 1:nPoles]
		curp2_i = [read(fid["Curp2$(j)"]) for j in 1:nPoles]
		close(fid)
	else
		curp1_i = [CurrentInit_PFree(hcen0,GammaT,chi1[i],kappa[i],beta,Gamma1) for i in 1:length(chi1)]
		curp2_i = [CurrentInit_PFree(hcen0,GammaT,chi2[i],kappa[i],beta,Gamma2) for i in 1:length(chi2)]
	end

	curM1		= CurrentMatrix(rho0,Gamma1,curp1_i)
	curM2		= CurrentMatrix(rho0,Gamma2,curp2_i)
	cur1_l		= [Current(curM1)]
 	cur2_l		= [Current(curM2)]
	scur1_l 	= [SpinCurrent(curM1)]
	scur2_l 	= [SpinCurrent(curM2)]

	
	for i in 1:numstep
		if resampling 
			if i % rsmpl_stp == 0
				A           = [SampleWignerA(beta,ws[i],vSampl) for i in 1:N]
				qin 		= [A[i][1] for i in 1:N]
				pin 		= [A[i][2] for i in 1:N]
			end
		end
		hcen_i  = FranssonHamiltonian(N,E,t,tso,a,c,num_laps,tph,tsoph,ws,qin)	
		chi1	= V[i]/2 .+ 1im .* xi ./beta 						# Overwrite previous chis by the new one WITH voltage drop
		chi2	= -V[i]/2 .+ 1im .* xi ./beta
	
		curp1_c, curp2_c, cur1_c, cur2_c, cur_c = RK_Currents(beta,hcen_i,Gamma1,Gamma2,GammaT,kappa,chi1,chi2,curp1_i,curp2_i,rho_i,dt)
		rho_c	= RK_VneumannL(rho_i,hcen_i,cur_c,dt)
		qn,pn   = SymplecticEulerStep(N,a,c,num_laps,tph,tsoph,ws,rho_c,qin,pin,dt) 
		# BUFFER
		rho_i   = copy(rho_c)
		curp1_i = copy(curp1_c)
		curp2_i = copy(curp2_c)
		qin     = qn
		pin     = pn
		pop_l,tim,cur1_l,cur2_l,scur1_l,scur2_l,q_l,p_l = SaveOutput(cur1_c,cur2_c,pop_l,rho_i,qn,pn,tim,i,dt,cur1_l,cur2_l,scur1_l,scur2_l,q_l,p_l,path)
		if i == numstep 
			fil = h5open("$(path)_final_Quantities","w")
			fil["Rho"] = rho_i
			for j in 1:nPoles
				fil["Curp1$(j)"] = curp1_i[j]
				fil["Curp2$(j)"] = curp2_i[j]
			end
			fil["q"] = qin
			fil["p"] = pin
			close(fil)	 
		end
	end
end

function PropagateSingle(E   		::Array,					 
						 ws  		::Float64,
						 gph 		::Number,
						 dt	 		::Number,
						 numstep	::Int,
						 tL			::Array,
						 tR			::Array,
						 beta		::Float64,
						 nPoles		::Int,
						 V          ::Number,
						 path		::String,
						 resampling ::Bool,
						 rhoIn      ::String,
						 phIn	    ::String,
						 vSampl     ::String,
						 pathIn     = "",
						 popIn 		= zeros(1),
						 qnIn  		= zeros(1),
						 pnIn  		= zeros(1),
						 rsmpl_stp	= 0)



	if phIn == "TDSCF"
		q0 			= qnIn
		p0 			= pnIn
	elseif phIn == "Wigner"
		A           = SampleWignerA(beta,ws,vSampl)
		q0 			= A[1]
		p0 			= A[2]
	elseif phIn == "fromFile"
		fid     = h5open(pathIn,"r")
		q0      = read(fid["q"])
		p0      = read(fid["p"])
	end

	
	qin         = q0
	pin         = p0
	q_l 		= [q0]
	p_l 		= [p0]
	
	hcen_i      = SingleLevel(E,ws,gph,q0)
	Gamma1,Gamma2 = GammaMatrices1d(1,tL,tR)
	GammaT		= Gamma1+Gamma2
	tim   		= [0.0]
	kappa,xi	= PadeFermi(nPoles)
	chi1		= V/2 .+ 1im .* xi ./beta    		# At t<t_0, voltage drop is zero.
	chi2		= -V/2 .+ 1im .* xi ./beta
	chiLP 	    = [chi1,chi2]

	if rhoIn == "pure"
		rho0    = Diagonal(float(popIn))
	elseif rhoIn == "equilibrium"
		rho0	    = DensityFExOpen(1,1,hcen_i,kappa,chiLP,GammaT,beta)*(1+0*1im)
	elseif rhoIn == "fromFile"
		rho0    = read(fid["Rho"])
	end
	
	pop_l       = [Populations(rho0)]
	rho_i       = rho0
	
	hcen0 		= hcen_i
	
	if rhoIn == "fromFile"
		curp1_i = [read(fid["Curp1$(j)"]) for j in 1:nPoles]
		curp2_i = [read(fid["Curp2$(j)"]) for j in 1:nPoles]
		close(fid)
	else
		curp1_i = [CurrentInit_PFree(hcen0,GammaT,chi1[i],kappa[i],beta,Gamma1) for i in 1:length(chi1)]
		curp2_i = [CurrentInit_PFree(hcen0,GammaT,chi2[i],kappa[i],beta,Gamma2) for i in 1:length(chi2)]
	end

	curM1		= CurrentMatrix(rho0,Gamma1,curp1_i)
	curM2		= CurrentMatrix(rho0,Gamma2,curp2_i)
	cur1_l		= [Current(curM1)]
	cur2_l		= [Current(curM2)]
	scur1_l 	= [SpinCurrent(curM1)]
	scur2_l 	= [SpinCurrent(curM2)]
	
	
	for i in 1:numstep
		if resampling 
			if i % rsmpl_stp == 0
				A           = SampleWignerA(beta,ws,vSampl) 
				qin 		= A[1]
				pin 		= A[2]
			end
		end
		hcen_i  = SingleLevel(E,ws,gph,qin)	
		chi1	= V/2 .+ 1im .* xi ./beta 						# Overwrite previous chis by the new one WITH voltage drop
		chi2	= -V/2 .+ 1im .* xi ./beta
	
		curp1_c, curp2_c, cur1_c, cur2_c, cur_c = RK_Currents(beta,hcen_i,Gamma1,Gamma2,GammaT,kappa,chi1,chi2,curp1_i,curp2_i,rho_i,dt)
		rho_c	= RK_VneumannL(rho_i,hcen_i,cur_c,dt)
		qn,pn   = SymplecticEulerStepSingle(ws,gph,rho_c,qin,pin,dt) 
		# BUFFER
		rho_i   = copy(rho_c)
		curp1_i = copy(curp1_c)
		curp2_i = copy(curp2_c)
		qin     = qn
		pin     = pn
		pop_l,tim,cur1_l,cur2_l,scur1_l,scur2_l,q_l,p_l = SaveOutput(cur1_c,cur2_c,pop_l,rho_i,qn,pn,tim,i,dt,cur1_l,cur2_l,scur1_l,scur2_l,q_l,p_l,path)
		if i == numstep 
			fil = h5open("$(path)_final_Quantities","w")
			fil["Rho"] = rho_i
			for j in 1:nPoles
				fil["Curp1$(j)"] = curp1_i[j]
				fil["Curp2$(j)"] = curp2_i[j]
			end
			fil["q"] = qin
			fil["p"] = pin
			close(fil)	 
		end
	end

end

end
