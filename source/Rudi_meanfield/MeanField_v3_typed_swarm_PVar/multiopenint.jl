include("global_constants.jl")
include("hamiltonian.jl")
include("iosave.jl")
include("observables.jl")
include("electronic_integration.jl")
include("classical_phonons.jl")
include("structures.jl")
using LinearAlgebra
using HDF5

function Propagate(
		  X_l		::Array,
		  P_l		::Array,
		  HamPar 	::Any,  
		  dt		::Number,
		  numstep	::Int,
		  tL		::Array{Float64,1},
		  tR		::Array{Float64,1},
		  beta		::Float64,
		  nPoles	::Int,
		  V         	::Float64,
		  path		::String,
		  rhoIn     	::String,
		  num_avg	::Int64,
		  elIntegrator 	::Function,
		  clIntegrator  ::Function,
		  popIn 	= zeros(1),
		  pathIn    	= "")


#	X0,P0 		= SampleWigner_Multi(beta,HamPar.ws,HamPar.N,num_avg)
#	X_l,P_l 	= zeros(num_avg,HamPar.N),zeros(num_avg,HamPar.N) 	
#	X_l[:,:] 	= X0
#	P_l[:,:]	= P0

#	HamPar.X = 1/num_avg * sum(X0,dims=1)	
#	HamPar.P = 1/num_avg * sum(P0,dims=1)

	X_avg_l 	= zeros(HamPar.N,LengthData)
	P_avg_l 	= zeros(HamPar.N,LengthData)

	hcen_i      	= HamPar.Hamiltonian(HamPar)
	Gamma1,Gamma2 	= GammaMatrices1d(HamPar.N,tL,tR)
	GammaT		= Gamma1+Gamma2
	tim	   	= zeros(LengthData)
	kappa,xi	= PadeFermi(nPoles)
	chi1		= 0.0 .+ 1im .* xi ./beta    		# At t<t_0, voltage drop is zero.
	chi2		= 0.0 .+ 1im .* xi ./beta
	chiLP 	    	= [chi1,chi2]
	if rhoIn == "pure"
		rho0    = zeros(ComplexF64,size(hcen_i)) .+ Diagonal(float(popIn))
	elseif rhoIn == "equilibrium"
		rho0	= DensityFExOpen(HamPar.N,1,hcen_i,kappa,chiLP,GammaT,beta)*(1+0*1im)
	elseif rhoIn == "closedEquilibrium"
		rho0    = DensityFEx(HamPar.N,1,hcen_i,HamPar.N,0)
	elseif rhoIn == "fromFile"
		rho0    = read(fid["Rho"])
	end
	
	conv_l   = zeros(3,LengthData)
	pop_l    = zeros(2*HamPar.N,LengthData)
	pop_es_l = zeros(2*HamPar.N,LengthData)
	rho_i	 = rho0
    
	
	if rhoIn == "fromFile"
		curp1_i = [read(fid["Curp1$(j)"]) for j in 1:nPoles]
		curp2_i = [read(fid["Curp2$(j)"]) for j in 1:nPoles]
		close(fid)
	else
		shx,shy = size(hcen_i)
		curp1_i = Array{ComplexF64}(undef,shx,shy,nPoles)
		curp2_i = Array{ComplexF64}(undef,shx,shy,nPoles)
		for i in 1:nPoles
			curp1_i[:,:,i] = CurrentInit_PFree(hcen_i,GammaT,chi1[i],kappa[i],beta,Gamma1)
			curp2_i[:,:,i] = CurrentInit_PFree(hcen_i,GammaT,chi2[i],kappa[i],beta,Gamma2)
		end
	end

	cur1_i		= CurrentMatrix(rho0,Gamma1,curp1_i)
	cur2_i		= CurrentMatrix(rho0,Gamma2,curp2_i)
	cur_i 		= cur1_i + cur2_i 
	cur1_l		= zeros(LengthData)
 	cur2_l		= zeros(LengthData)
	scur1_l 	= zeros(3,LengthData)
	scur2_l 	= zeros(3,LengthData)
	mom_el      	= zeros(2*(HamPar.N-1),LengthData)
	loc_cur_l   	= zeros(HamPar.N-1,LengthData)

	chi1		= V/2 .+ 1im .* xi ./beta 						# Overwrite previous chis by the new one WITH voltage drop
	chi2		= -V/2 .+ 1im .* xi ./beta
	vind 		= 1	
	for i in 1:numstep
				
		Xn,Pn = clIntegrator(HamPar,rho_i,dt,X_l,P_l,num_avg)
		X_l   = Xn
		P_l   = Pn
		HamPar.X = 1/num_avg * sum(Xn,dims=1)
		HamPar.P = 1/num_avg * sum(Pn,dims=1)
		
		cur1_i,cur2_i = elIntegrator(shx,shy,nPoles,beta,hcen_i,Gamma1,Gamma2,GammaT,kappa,chi1,chi2,curp1_i,curp2_i,cur_i,rho_i,dt)
		hcen_i  = HamPar.Hamiltonian(HamPar)	
		conv_l,pop_l,pop_es_l,tim,cur1_l,cur2_l,scur1_l,scur2_l,X_avg_l,P_avg_l,mom_el,loc_cur_l,vind = SaveOutput(HamPar,cur1_i,cur2_i,conv_l,pop_l,pop_es_l,rho_i,hcen_i,tim,i,dt,cur1_l,cur2_l,scur1_l,scur2_l,loc_cur_l,X_avg_l,P_avg_l,mom_el,path,vind)
		if i == numstep 
			fil = h5open("$(path)_final_Quantities","w")
			fil["Rho"] = rho_i
			for j in 1:nPoles
				fil["Curp1$(j)"] = curp1_i[j]
				fil["Curp2$(j)"] = curp2_i[j]
			end
			fil["X"] = float.(HamPar.X)
			fil["P"] = float.(HamPar.P)
			close(fil)	 
		end
	end
end

