module IOsave
include("electrons.jl")
include("global_constants.jl")
using .GlobalConstants
using .Electrons
using LinearAlgebra
using Dates
using DelimitedFiles
export SaveParams,SaveOutput,Average

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

function SaveOutput(cur1_c,cur2_c,pop_l,pop_es_l,rho_i,hcen_i,qn,pn,tim,i,dt,cur1_l,cur2_l,loc_cur_l,q_l,p_l,mom_el,weight,path)

	pop = Populations(rho_i)
	pop_es = PopulationsES(rho_i,hcen_i)
	mom = MomentumElectrons(rho_i)
	cur1 = Current(cur1_c)
	cur2 = Current(cur2_c)
	loccur = LocalCurrent(rho_i,hcen_i)
	if i == 1
		curfile = open("$(path)_currents.dat","w")
		popfile = open("$(path)_electrons.dat","w")
		popesfile = open("$(path)_electrons_es.dat","w")
		phfile  = open("$(path)_phonons.dat","w")
		weightf = open("$(path)_weight.dat","w")
		momentumf= open("$(path)_momentum_el.dat","w")
		loccurf = open("$(path)_loc_cur.dat","w")

		write(curfile,"#(1)t (2)IL (3)IR\n")	
		write(popfile,"#(1)t ")
		write(popesfile,"#(1)t ")
		for k in 1:length(qn)
			write(popfile,"($(k+1))p$(k) ")
		end
		for k in 1:length(pop_es)
                        write(popesfile,"($(k+1))S$(k) ")
                end	
		
		write(popfile,"\n")
		write(popesfile,"\n")
		write(phfile,"#(1)t ")
		k = 2
		for i in 1:length(qn)
			write(phfile,"($(k))q$(i) ($(k+1))p$(i) ")
			k += 2
		end
		write(phfile,"\n")
		write(weightf,"$(weight)\n")

		write(momentumf,"#(1)t ")
		k = 2
		for i in 1:length(qn)-1
			write(momentumf,"($(i+1))p$(i) ")
		end	
		write(momentumf,"\n")

		write(loccurf,"#(1)t")
		for i in 1:length(qn)-1
			write(loccurf,"($(i+1))cur($(i),$(i+1)) ")
		end
		write(loccurf,"\n")

		close(curfile)
		close(popfile)
		close(popesfile)
		close(phfile)
		close(weightf)
		close(momentumf)
		close(loccurf)
	end
	
	if i % GlobalConstants.AppendData == 0
		push!(pop_l,pop)
		push!(pop_es_l,pop_es)
		push!(cur1_l,cur1)
		push!(cur2_l,cur2)
		push!(q_l, qn)
		push!(p_l, pn)
		push!(mom_el, mom)
		push!(loc_cur_l,loccur)
		push!(tim,abs(i*dt))
	end

	if i % GlobalConstants.SaveData == 0
		curfile = open("$(path)_currents.dat","a")
		popfile = open("$(path)_electrons.dat","a")
		popesfile = open("$(path)_electrons_es.dat","a")
		phfile  = open("$(path)_phonons.dat","a")
		momfile = open("$(path)_momentum_el.dat","a")
		loccurf = open("$(path)_loc_cur.dat","a")
		for j in 1:length(pop_l)
			write(popfile,"$(tim[j]) ")
			write(popesfile,"$(tim[j]) ")
			write(curfile,"$(tim[j]) $(cur1_l[j]) $(cur2_l[j]) \n") 
			write(phfile, "$(tim[j]) ")
			write(momfile, "$(tim[j]) ")
			write(loccurf, "$(tim[j]) ")

				
			for k in 1:length(pop_l[j])
				write(popfile,"$(pop_l[j][k]) ")
				write(popesfile,"$(pop_es_l[j][k]) ")
			end
			for k in 1:length(q_l[j])
				write(phfile, "$(q_l[j][k]) $(p_l[j][k]) ")
			end
			for k in 1:length(mom_el[j])
				write(momfile, "$(mom_el[j][k]) ")
			end
			for k in 1:length(loc_cur_l[j])
				write(loccurf, "$(loc_cur_l[j][k]) ")
			end

			write(popfile, "\n")
			write(popesfile, "\n")
			write(phfile, "\n")
			write(momfile, "\n")
			write(loccurf, "\n")
		end
		close(curfile)
		close(popfile)		
		close(popesfile)	
		close(phfile)
		close(momfile)
		close(loccurf)
		pop_l = []
		pop_es_l = []
		cur1_l = []
		cur2_l = []
		q_l = []
		p_l = []
		mom_el = []
		loc_cur_l = []
		tim = []

	end


	return pop_l,pop_es_l,tim,cur1_l,cur2_l,q_l,p_l,mom_el,loc_cur_l
end

function Average(ind,path,num_avg)
	pathr = "$(path)_$(ind)"
	observables = ["currents","electrons","electrons_es","phonons","momentum_el","loc_cur"]
	rm("$(pathr)_final_Quantities")
	if ind == 1
		weight  = readdlm("$(pathr)_weight.dat")[1]	
		weightTf = open("$(path)_weightTot.dat","w")
		write(weightTf,"$(weight)")
		close(weightTf)
		rm("$(pathr)_weight.dat")
		
		for obs in observables
			file_avg = readdlm("$(pathr)_$(obs).dat")
			file_avg[2:end,2:end] = file_avg[2:end,2:end] * weight
			writedlm("$(path)_$(obs)_avg.dat",file_avg) 
			rm("$(pathr)_$(obs).dat")
		end
	elseif ind < num_avg
		weight  = readdlm("$(pathr)_weight.dat")[1]
		weightT = readdlm("$(path)_weightTot.dat")[1]	
		weightT += weight

		writedlm("$(path)_weightTot.dat",weightT) 
		rm("$(pathr)_weight.dat")
					
		for obs in observables
			file_avg = readdlm("$(path)_$(obs)_avg.dat")
			file_raw = readdlm("$(pathr)_$(obs).dat")

			F_avg    = copy(file_avg)
			F_avg[2:end,2:end] += file_raw[2:end,2:end] * weight
			rm("$(path)_$(obs)_avg.dat")
			writedlm("$(path)_$(obs)_avg.dat",F_avg) 
			rm("$(pathr)_$(obs).dat")
		end
	elseif ind == num_avg
		weight  = readdlm("$(pathr)_weight.dat")[1]	
		weightT = readdlm("$(path)_weightTot.dat")[1]	
		weightT += weight

		writedlm("$(path)_weightTot.dat",weightT) 
		rm("$(pathr)_weight.dat")
		for obs in observables
			file_avg = readdlm("$(path)_$(obs)_avg.dat")
			file_raw = readdlm("$(pathr)_$(obs).dat")

			F_avg    = copy(file_avg)
			F_avg[2:end,2:end] += file_raw[2:end,2:end] * weight 
			F_avg[2:end,2:end] /= weightT
			rm("$(path)_$(obs)_avg.dat")
			writedlm("$(path)_$(obs)_avg.dat",F_avg) 
			rm("$(pathr)_$(obs).dat")
		end
		
	end
end

end
