include("observables.jl")
include("global_constants.jl")
using LinearAlgebra
using Dates
using DelimitedFiles

function SaveParams(params,path,purpose)
	file = open("$(path)","w")
	write(file,"# Parameters of run on $(Dates.now()) \n")
	write(file,"# The purpose of this run is the following: $(purpose)\n")
	write(file,"# --------------------PARAMS----------------------------\n")
	for pars in params
		write(file,"$(pars[1]) = $(pars[2])\n")
	end
	write(file,"# ----------------THIS IS THE END.----------------------\n")	
	close(file)
end

function SaveOutput(HamPar::Any,
		    cur1_c::Array{ComplexF64},
		    cur2_c::Array{ComplexF64},
		    conv_l::Array{Float64},
		    pop_l::Array{Float64},
		    pop_es_l::Array{Float64},
		    rho_i::Array{ComplexF64},
		    hcen_i::Array{ComplexF64},
		    tim::Array{Float64},
		    i::Int64,
		    dt::Float64,
		    cur1_l::Array{Float64},
		    cur2_l::Array{Float64},
		    scur1_l::Array{Float64},
		    scur2_l::Array{Float64},
		    loc_cur_l::Array{Float64},
		    q_l::Array{Float64},
		    p_l::Array{Float64},
		    mom_el::Array{Float64},
		    path::String,
		    vind::Int64)::Tuple{Array{Float64},Array{Float64},Array{Float64},Array{Float64},Array{Float64},Array{Float64},Array{Float64},Array{Float64},Array{Float64},Array{Float64},Array{Float64},Array{Float64},Int64}

	if i == 1
		convfile= open("$(path)_convergence_info.dat","w")
		curfile = open("$(path)_currents.dat","w")
		popfile = open("$(path)_electrons.dat","w")
		popesfile = open("$(path)_electrons_es.dat","w")
		phfile  = open("$(path)_phonons.dat","w")
		momentumf= open("$(path)_momentum_el.dat","w")
		loccurf = open("$(path)_loc_cur.dat","w")
		
		write(convfile,"#(1)t (2)trace (3)sys_energy (4)dif\n")
		write(curfile,"#(1)t (2)IL (3)IR (4)sILx (5)sILy (6)sILz (7)sIRx (8)sIRy (9)sIRz\n")	
		write(popfile,"#(1)t ")
		write(popesfile,"#(1)t ")
		for k in 1:HamPar.N
			write(popfile,"($(k+1))p$(k)u ($(k+1))p$(k)d ")
		end
		for k in 1:Int(2*HamPar.N)
                        write(popesfile,"($(k+1))S$(k) ")
                end	
		
		write(popfile,"\n")
		write(popesfile,"\n")
		write(phfile,"#(1)t ")
		k = 2
		for i in 1:HamPar.N
			write(phfile,"($(k))q$(i) ($(k+1))p$(i) ")
			k += 2
		end
		write(phfile,"\n")

		write(momentumf,"#(1)t ")
		k = 2
		for i in 1:HamPar.N-1
			write(momentumf,"($(k))p$(i)u ($(k+1))p$(i)d ")
			k += 2
		end	
		write(momentumf,"\n")

		write(loccurf,"#(1)t")
		for i in 1:HamPar.N-1
			write(loccurf,"($(i+1))cur($(i),$(i+1)) ")
		end
		write(loccurf,"\n")

		close(convfile)
		close(curfile)
		close(popfile)
		close(popesfile)
		close(phfile)
		close(momentumf)
		close(loccurf)
	end
	
	if i % AppendData == 0

		conv_l[:,vind] = ConvergenceInfo(HamPar,rho_i,hcen_i)
		pop_l[:,vind] = Populations(rho_i)
		pop_es_l[:,vind]= PopulationsES(rho_i,hcen_i)
		cur1_l[vind] = Current(cur1_c)
		cur2_l[vind] = Current(cur2_c)
		scur1_l[:,vind] = SpinCurrent(cur1_c)
		scur2_l[:,vind] = SpinCurrent(cur2_c)
		q_l[:,vind] = HamPar.X
		p_l[:,vind] = HamPar.P
		mom_el[:,vind] = MomentumElectrons(rho_i)
		loc_cur_l[:,vind] = LocalCurrent(rho_i,hcen_i)
		tim[vind] = abs(i*dt)
		vind += 1
	end

	if i % SaveData == 0
		convfile = open("$(path)_convergence_info.dat","a")		
		curfile = open("$(path)_currents.dat","a")
		popfile = open("$(path)_electrons.dat","a")
		popesfile = open("$(path)_electrons_es.dat","a")
		phfile  = open("$(path)_phonons.dat","a")
		momfile = open("$(path)_momentum_el.dat","a")
		loccurf = open("$(path)_loc_cur.dat","a")
		for j in 1:length(tim)
			write(convfile,"$(tim[j]) $(real(conv_l[1,j])) $(real(conv_l[2,j])) $(real(conv_l[3,j])) \n")
			write(popfile,"$(tim[j]) ")
			write(popesfile,"$(tim[j]) ")
			write(curfile,"$(tim[j]) $(cur1_l[j]) $(cur2_l[j]) $(scur1_l[1,j]) $(scur1_l[2,j]) $(scur1_l[3,j]) $(scur2_l[1,j]) $(scur2_l[2,j]) $(scur2_l[3,j]) \n")
			write(phfile, "$(tim[j]) ")
			write(momfile, "$(tim[j]) ")
			write(loccurf, "$(tim[j]) ")

				
			for k in 1:length(pop_l[:,j])
				write(popfile,"$(pop_l[k,j]) ")
				write(popesfile,"$(pop_es_l[k,j]) ")
			end
			for k in 1:length(q_l[:,j])
				write(phfile, "$(q_l[k,j]) $(p_l[k,j]) ")
			end
			for k in 1:length(mom_el[:,j])
				write(momfile, "$(mom_el[k,j]) ")
			end
			for k in 1:length(loc_cur_l[:,j])
				write(loccurf, "$(loc_cur_l[k,j]) ")
			end

			write(popfile, "\n")
			write(popesfile, "\n")
			write(phfile, "\n")
			write(momfile, "\n")
			write(loccurf, "\n")
		end
		close(convfile)
		close(curfile)
		close(popfile)		
		close(popesfile)	
		close(phfile)
		close(momfile)
		close(loccurf)
		vind = 1
	end

	return conv_l,pop_l,pop_es_l,tim,cur1_l,cur2_l,scur1_l,scur2_l,q_l,p_l,mom_el,loc_cur_l,vind
end

function Average(ind,path,num_avg)
	pathr = "$(path)_$(ind)"
	observables = ["currents","electrons","electrons_es","phonons","momentum_el","loc_cur","convergence_info"]
	rm("$(pathr)_final_Quantities")
	if ind == 1
		for obs in observables
			file_avg = readdlm("$(pathr)_$(obs).dat")
			writedlm("$(path)_$(obs)_avg.dat",file_avg) 
			rm("$(pathr)_$(obs).dat")
		end
	elseif ind <= num_avg
	
		for obs in observables
			file_avg = readdlm("$(path)_$(obs)_avg.dat")
			file_raw = readdlm("$(pathr)_$(obs).dat")

			F_avg    = copy(file_avg)
			F_avg[2:end,2:end] += file_raw[2:end,2:end]
			F_avg[2:end,2:end] = F_avg[2:end,2:end]
			rm("$(path)_$(obs)_avg.dat")
			writedlm("$(path)_$(obs)_avg.dat",F_avg) 
			rm("$(pathr)_$(obs).dat")
		end
	elseif ind == num_avg
		for obs in observables
                        file_avg = readdlm("$(path)_$(obs)_avg.dat")
                        file_raw = readdlm("$(pathr)_$(obs).dat")

                        F_avg    = copy(file_avg)
                        F_avg[2:end,2:end] .+= file_raw[2:end,2:end]
                        F_avg[2:end,2:end] ./= F_avg[2:end,2:end]
                        rm("$(path)_$(obs)_avg.dat")
                        writedlm("$(path)_$(obs)_avg.dat",F_avg)
                        rm("$(pathr)_$(obs).dat")
                end
	end
end

