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
		  HamPar 	::Any,  
		  dt		::Number,
		  numstep	::Int,
		  tL		::Array{Float64,1},
		  tR		::Array{Float64,1},
		  beta		::Float64,
		  nPoles	::Int,
		  V         	::Float64,
		  path		::String,
		  elIntegrator 	::Function,
		  clIntegrator  ::Function,
		  initProp	::Bool,
		  rhoIn     	::String,
		  popIn 	::Array{Float64,1},
		  pathIn    	::String,
		  propStp	::Int)

	if initProp 
		X_l 	    	= zeros(HamPar.N,LengthData)
		P_l 	    	= zeros(HamPar.N,LengthData)

		hcen_i      	= HamPar.Hamiltonian(HamPar)*(1+0.0*1im)
		tim	   	= zeros(LengthData)
		kappa,xi	= PadeFermi(nPoles)
		chi1		= 0.0 .+ 1im .* xi ./beta    		# At t<t_0, voltage drop is zero.
		chi2		= 0.0 .+ 1im .* xi ./beta
		chiLP 	    	= [chi1,chi2]


		conv_l   = zeros(3,LengthData)
		pop_l    = zeros(2*HamPar.N,LengthData)
		pop_es_l = zeros(2*HamPar.N,LengthData)
    
		shx,shy = size(hcen_i)
		curp1_i = Array{ComplexF64}(undef,shx,shy,nPoles)
		curp2_i = Array{ComplexF64}(undef,shx,shy,nPoles)
		Gamma1,Gamma2 	= GammaMatrices1d(HamPar.N,tL,tR)
		GammaT		= Gamma1+Gamma2
		if rhoIn == "pure"
			rho0    = zeros(ComplexF64,size(hcen_i)) .+ Diagonal(float(popIn))
		elseif rhoIn == "equilibrium"
			rho0	= DensityFExOpen(HamPar.N,1,hcen_i,kappa,chiLP,GammaT,beta)*(1+0*1im)
		elseif rhoIn == "closedEquilibrium"
			rho0    = DensityFEx(HamPar.N,1,hcen_i,HamPar.N,0)
		elseif rhoIn == "fromFile"
			fid 	= h5open(pathIn)
			rho0    = read(fid["Rho"])
		end
	
		if rhoIn == "fromFile"
			fid     = h5open(pathIn)
			for i in 1:nPoles	
				curp1_i[:,:,i] .= read(fid["Curp1$(i)"]) 
				curp2_i[:,:,i] .= read(fid["Curp2$(i)"]) 
			end
			close(fid)
		else
			for i in 1:nPoles
				curp1_i[:,:,i] = CurrentInit_PFree(hcen_i,GammaT,chi1[i],kappa[i],beta,Gamma1)
				curp2_i[:,:,i] = CurrentInit_PFree(hcen_i,GammaT,chi2[i],kappa[i],beta,Gamma2)
			end
		end
		
		rho_i	 	= rho0
		cur1_i		= CurrentMatrix(rho0,Gamma1,curp1_i)
		cur2_i		= CurrentMatrix(rho0,Gamma2,curp2_i)
		cur_i 		= cur1_i + cur2_i 
		cur1_l		= zeros(LengthData)
		cur2_l		= zeros(LengthData)
		scur1_l 	= zeros(3,LengthData)
		scur2_l 	= zeros(3,LengthData)

		tInit = 1
		tEnd = numstep	
		vind 		= 1
	else
		kappa,xi        = PadeFermi(nPoles)	
		fid     = h5open(pathIn)
		HamPar.X = read(fid["X"])
		HamPar.P = read(fid["P"])
		hcen_i = HamPar.Hamiltonian(HamPar)*(1+0.0*1im)	
		rho0 = read(fid["Rho"])
		shx,shy = size(hcen_i)
                curp1_i = Array{ComplexF64}(undef,shx,shy,nPoles)
                curp2_i = Array{ComplexF64}(undef,shx,shy,nPoles)
		for i in 1:nPoles
			curp1_i[:,:,i] .= read(fid["Curp1$(i)"])
			curp2_i[:,:,i] .= read(fid["Curp2$(i)"])
		end
		X_l = read(fid["X_l"])	
		P_l = read(fid["P_l"])
		tim = read(fid["tim"])
		conv_l   = read(fid["conv_l"])
                pop_l    = read(fid["pop_l"])
                pop_es_l = read(fid["pop_es_l"])
		rho_i    = rho0
		Gamma1,Gamma2   = GammaMatrices1d(HamPar.N,tL,tR)
                GammaT          = Gamma1+Gamma2
		cur1_i   = CurrentMatrix(rho0,Gamma1,curp1_i)
		cur2_i   = CurrentMatrix(rho0,Gamma2,curp2_i)
		cur_i    = cur1_i + cur2_i
		cur1_l   = read(fid["cur1_l"])
		cur2_l   = read(fid["cur2_l"])
		scur1_l  = read(fid["scur1_l"])
                scur2_l  = read(fid["scur2_l"])

		tInit    = propStp*numstep + 1
		tEnd     = (tInit-1) + numstep
		vind     = read(fid["vind"])
	end
	
#	mom_el      	= zeros(2*(HamPar.N-1),LengthData)
#	loc_cur_l   	= zeros(HamPar.N-1,LengthData)

	chi1		= V/2 .+ 1im .* xi ./beta 						# Overwrite previous chis by the new one WITH voltage drop
	chi2		= -V/2 .+ 1im .* xi ./beta
	
	
	for i in tInit:tEnd
				
		HamPar.X,HamPar.P = clIntegrator(HamPar,rho_i,dt)
		cur1_i,cur2_i = elIntegrator(shx,shy,nPoles,beta,hcen_i,Gamma1,Gamma2,GammaT,kappa,chi1,chi2,curp1_i,curp2_i,cur_i,rho_i,dt)
		hcen_i  = HamPar.Hamiltonian(HamPar)*(1+0.0*1im)	
		conv_l,pop_l,pop_es_l,tim,cur1_l,cur2_l,scur1_l,scur2_l,X_l,P_l,vind = SaveOutput(HamPar,cur1_i,cur2_i,conv_l,pop_l,pop_es_l,rho_i,hcen_i,tim,i,dt,cur1_l,cur2_l,scur1_l,scur2_l,X_l,P_l,path,vind)
		if i == tEnd 
			fil = h5open("$(path)_final_Quantities","w")
			fil["Rho"] = rho_i
			for j in 1:nPoles
				fil["Curp1$(j)"] = curp1_i[j]
				fil["Curp2$(j)"] = curp2_i[j]
			end
			fil["X"] = float.(HamPar.X)
			fil["P"] = float.(HamPar.P)
			fil["X_l"] = X_l
                        fil["P_l"] = P_l
			fil["tim"] = tim
			fil["conv_l"] = conv_l
			fil["pop_l"] = pop_l
			fil["pop_es_l"] = pop_es_l
			fil["cur1_l"] = cur1_l
			fil["cur2_l"] = cur2_l
			fil["scur1_l"] = scur1_l
                        fil["scur2_l"] = scur2_l
			fil["vind"] = vind
			close(fil)	 
		end
	end
end

