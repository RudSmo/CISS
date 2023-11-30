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
		  HamPar::Any,  
		  dt		::Number,
		  numstep	::Int,
		  tL		::Array{Float64,1},
		  tR		::Array{Float64,1},
		  beta		::Float64,
		  nPoles	::Int,
		  V         	::Float64,
		  path		::String,
		  rhoIn     	::String,
		  elIntegrator 	::Function,
		  clIntegrator  ::Function,
		  initProp      ::Bool,
		  Lauf			= 0, # Laufindex der die Propagationen zählt von 0 bis n 
		  popIn 	= zeros(1),
		  pathIn    	= ""
		 
		  )


#                       ???????? wo wird LengthData definiert???????????? -> global constants
	
	#loc_cur_l   = zeros(HamPar.N-1,LengthData)

	if initProp
		X_l 	 = zeros(HamPar.N,LengthData) 
		P_l 	 = zeros(HamPar.N,LengthData)
		conv_l   = zeros(3,LengthData)
		pop_l    = zeros(2*HamPar.N,LengthData)
		pop_es_l = zeros(2*HamPar.N,LengthData)

		hcen_i      	= HamPar.Hamiltonian(HamPar)*(1+0.0*1im)
		Gamma1,Gamma2 	= GammaMatrices1d(HamPar.N,tL,tR)
		GammaT		= Gamma1+Gamma2
		tim	   	= zeros(LengthData)
		kappa,xi	= PadeFermi(nPoles)
		chi1		= 0.0 .+ 1im .* xi ./beta    		# At t<t_0, voltage drop is zero.
		chi2		= 0.0 .+ 1im .* xi ./beta
		chiLP 	   	= [chi1,chi2]

		shx,shy = size(hcen_i)
		curp1_i = Array{ComplexF64}(undef,shx,shy,nPoles)
		curp2_i = Array{ComplexF64}(undef,shx,shy,nPoles)

		cur1_l		= zeros(LengthData)
		cur2_l		= zeros(LengthData)
		scur1_l 	= zeros(3,LengthData)
		scur2_l 	= zeros(3,LengthData)
		if rhoIn == "pure" # in input_Sophia kann rhoIn = pure sein oder closedequilib ??
			rho0    = zeros(ComplexF64,size(hcen_i)) .+ Diagonal(float(popIn))
		elseif rhoIn == "equilibrium"
			rho0	= DensityFExOpen(HamPar.N,1,hcen_i,kappa,chiLP,GammaT,beta)*(1+0*1im)
		elseif rhoIn == "closedEquilibrium"
			rho0    = DensityFEx(HamPar.N,1,hcen_i,HamPar.N,0)
		elseif rhoIn == "fromFile"
			#pathin istauch data_final_quantities = path
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
		cur1_i		= CurrentMatrix(rho0,Gamma1,curp1_i)
		cur2_i		= CurrentMatrix(rho0,Gamma2,curp2_i)
		cur_i 		= cur1_i + cur2_i 
		chi1		= V/2 .+ 1im .* xi ./beta 						# Overwrite previous chis by the new one WITH voltage drop
		chi2		= -V/2 .+ 1im .* xi ./beta
		#was ist vind nochmal für ein Laufindex?
		vind 		= 1	

	else # 
		last_fil = h5open("$(path)_final_Quantities","r") 
		HamPar.X = read(last_fil["X"])
		HamPar.P = read(last_fil["P"])

		hcen_i      	= HamPar.Hamiltonian(HamPar)*(1+0.0*1im)
		Gamma1,Gamma2 	= GammaMatrices1d(HamPar.N,tL,tR)
		GammaT		= Gamma1+Gamma2
		
		kappa,xi	= PadeFermi(nPoles)
		shx,shy = size(hcen_i)
		curp1_i = Array{ComplexF64}(undef,shx,shy,nPoles)
		curp2_i = Array{ComplexF64}(undef,shx,shy,nPoles)

		# anderer pfad als lauf = 0
		# der letzte wert für rho, und curp.. soll als neuer anfangswert verwendet werden -> daten werden überschrieben
		rho0    = read(last_fil["Rho"]) 
		for i in 1:nPoles	
			curp1_i[:,:,i] .= read(last_fil["Curp1$(i)"]) 
			curp2_i[:,:,i] .= read(last_fil["Curp2$(i)"]) 
		end
		conv_l = read(last_fil["convergence"])
		pop_l = read(last_fil["populations"])
		pop_es_l = read(last_fil["populations_es"])
		# soll cur_i so. was sind dasfür ströme?
		cur1_i = CurrentMatrix(rho0, Gamma1, curp1_i)
		cur2_i = CurrentMatrix(rho0, Gamma2, curp2_i) #?=#CurrentMatrix(rho0, Gamma1, curp1_i)
		cur_i = cur1_i + cur2_i
		chi1 = V/2 .+ 1im .* xi ./ beta  # Overwrite previous chis by the new one WITH voltage drop
		chi2 = -V/2 .+ 1im .* xi ./ beta

		tim = read(last_fil["tim"])
		cur1_l = read(last_fil["currents1"])
		cur2_l = read(last_fil["currents2"])
		scur1_l = read(last_fil["spin_currents1"])
		scur2_l = read(last_fil["spin_currents2"])
		X_l = read(last_fil["xl"])
		P_l = read(last_fil["pl"])
		
		#was ist vind nochmal für ein Laufindex?
		#???????????????
		vind = read(last_fil["vind"])
		close(last_fil)
	end

	

	for i in 1+Lauf*numstep:(Lauf+1)*numstep
		
		HamPar.X,HamPar.P = clIntegrator(HamPar,rho0,dt)
		cur1_i,cur2_i = elIntegrator(shx,shy,nPoles,beta,hcen_i,Gamma1,Gamma2,GammaT,kappa,chi1,chi2,curp1_i,curp2_i,cur_i,rho0,dt)
		hcen_i  = HamPar.Hamiltonian(HamPar)*(1+0.0*1im)
		
		
		conv_l,pop_l,pop_es_l,tim,cur1_l,cur2_l,scur1_l,scur2_l,X_l,P_l,vind = SaveOutput(HamPar,cur1_i,cur2_i,conv_l,pop_l,pop_es_l,rho0,hcen_i,tim,i,dt,cur1_l,cur2_l,scur1_l,scur2_l,X_l,P_l,path,vind)#loc_cur_l,
		if i == (Lauf + 1)* numstep		
			fil = h5open("$(path)_final_Quantities","w") 
			fil["Rho"] = rho0
			for j in 1:nPoles
				fil["Curp1$(j)"] = curp1_i[j]
				fil["Curp2$(j)"] = curp2_i[j]
			end
			fil["X"] = float.(HamPar.X)
			fil["P"] = float.(HamPar.P)

			fil["convergence"] = conv_l
			fil["populations"]=  pop_l
			fil["populations_es"] =pop_es_l
			fil["tim"]=tim
			fil["currents1"] =cur1_l
			fil["currents2"]=cur2_l
			fil["spin_currents1"]=scur1_l
			fil["spin_currents2"]=scur2_l
			fil["xl"]=X_l
			fil["pl"]=P_l
			fil["vind"] =vind
			close(fil)
		end

		
	end
end

