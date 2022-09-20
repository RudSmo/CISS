module MultiOpenIntegrator
include("global_constants.jl")
include("hamiltonian.jl")
include("iosave.jl")
include("electrons.jl")
include("rungekutta.jl")
using .GlobalConstants
using .Hamiltonian
using .IOsave
using .Electrons
using .RungeKutta
using LinearAlgebra
using HDF5

export Propagate


function Propagate(N 		::Int,
		  hcen   	::Array,
		  dt		::Number,
		  numstep	::Int,
		  tL		::Array,
		  tR		::Array,
		  beta		::Float64,
		  nPoles	::Int,
		  muL           ::Float64,
		  muR           ::Float64,
		  path		::String,
	          rhoIn         ::String,
		  popIn         ::Array,
		  pathIn        = "")

	Gamma1,Gamma2   = GammaMatrices1d(N,tL,tR)
	GammaT		= Gamma1+Gamma2
	tim   		= [0.0]
	kappa,xi	= PadeFermi(nPoles)
	V0              = 0.0
	chi1		= V0/2 .+ 1im .* xi ./beta    		# At t<t_0, voltage drop is zero.
	chi2		= -V0/2 .+ 1im .* xi ./beta
	chiLP 	        = [chi1,chi2]

	if rhoIn == "pure"
		rho_i    = Diagonal(float(popIn))
	elseif rhoIn == "equilibriumOpen"
		rho_i	 = DensityFExOpen(N,1,hcen,kappa,chiLP,GammaT,beta)*(1+0*1im)
	elseif rhoIn == "equilibriumClosed"
		rho_i    = DensityFEx(N,1,hcen,N,0) 
	elseif rhoIn == "fromFile"
		fid      = h5open(pathIn,"r")
		rho_i    = read(fid["Rho"])
	end

	pop_l       	= [Populations(rho_i)]

	if rhoIn == "fromFile"
		curp1_i = [read(fid["Curp1$(j)"]) for j in 1:nPoles]
		curp2_i = [read(fid["Curp2$(j)"]) for j in 1:nPoles]
		close(fid)
	else
		curp1_i = [CurrentInit_PFree(hcen,GammaT,chi1[i],kappa[i],beta,Gamma1) for i in 1:length(chi1)]
		curp2_i = [CurrentInit_PFree(hcen,GammaT,chi2[i],kappa[i],beta,Gamma2) for i in 1:length(chi2)]
	end

	curM1		= CurrentMatrix(rho_i,Gamma1,curp1_i)
	curM2		= CurrentMatrix(rho_i,Gamma2,curp2_i)
	
	cur1_l		= [Current(curM1)]
	cur1u_l         = [CurrentSP(curM1)[1]]
	cur1d_l         = [CurrentSP(curM1)[2]]

 	cur2_l		= [Current(curM2)]
	cur2u_l         = [CurrentSP(curM2)[1]]
        cur2d_l         = [CurrentSP(curM2)[2]]

	scur1_l         = [SpinCurrent(curM1)]
	scur2_l         = [SpinCurrent(curM2)]
	
	for i in 1:numstep
		chi1	= muL .+ 1im .* xi ./beta 						# Overwrite previous chis by the new one WITH voltage drop
		chi2	= muR .+ 1im .* xi ./beta
	
		curp1_c, curp2_c, cur1_c, cur2_c, cur_c = RK_Currents(beta,hcen,Gamma1,Gamma2,GammaT,kappa,chi1,chi2,curp1_i,curp2_i,rho_i,dt)
		rho_c	= RK_VneumannL(rho_i,hcen,cur_c,dt)
		# BUFFER
		rho_i   = copy(rho_c)
		curp1_i = copy(curp1_c)
		curp2_i = copy(curp2_c)
#		pop_l,tim,cur1_l,cur2_l,scur1_l,scur2_l = SaveOutput(cur1_c,cur2_c,pop_l,rho_i,tim,i,dt,cur1_l,cur2_l,scur1_l,scur2_l,path)
		pop_l,tim,cur1_l,cur2_l,cur1u_l,cur1d_l,cur2u_l,cur2d_l,scur1_l,scur2_l = SaveOutput(cur1_c,cur2_c,pop_l,rho_i,tim,i,dt,cur1_l,cur2_l,cur1u_l,cur1d_l,cur2u_l,cur2d_l,scur1_l,scur2_l,path)
		if i == numstep 
			fil = h5open("$(path)_final_Quantities","w")
			fil["Rho"] = rho_i
			for j in 1:nPoles
				fil["Curp1$(j)"] = curp1_i[j]
				fil["Curp2$(j)"] = curp2_i[j]
			end
			close(fil)	 
		end
	end
end

end
