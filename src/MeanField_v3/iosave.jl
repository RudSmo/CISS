module IOsave
include("observables.jl")
include("global_constants.jl")
using .GlobalConstants
using .Observables
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

function SaveOutput(HamPar,cur1_c,cur2_c,conv_l,pop_l,pop_es_l,rho_i,hcen_i,qn,pn,tim,i,dt,cur1_l,cur2_l,scur1_l,scur2_l,loc_cur_l,q_l,p_l,mom_el,path)

	pop = Populations(rho_i)
	pop_es = PopulationsES(rho_i,hcen_i)
	mom = MomentumElectrons(rho_i)
	cur1 = Current(cur1_c)
	cur2 = Current(cur2_c)
	scur1 = SpinCurrent(cur1_c)
	scur2 = SpinCurrent(cur2_c)
	loccur = LocalCurrent(rho_i,hcen_i)
	conv = ConvergenceInfo(HamPar,rho_i,hcen_i) 
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
		for k in 1:length(qn)
			write(popfile,"($(k+1))p$(k)u ($(k+1))p$(k)d ")
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

		write(momentumf,"#(1)t ")
		k = 2
		for i in 1:length(qn)-1
			write(momentumf,"($(k))p$(i)u ($(k+1))p$(i)d ")
			k += 2
		end	
		write(momentumf,"\n")

		write(loccurf,"#(1)t")
		for i in 1:length(qn)-1
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
	
	if i % GlobalConstants.AppendData == 0
		push!(conv_l,conv)
		push!(pop_l,pop)
		push!(pop_es_l,pop_es)
		push!(cur1_l,cur1)
		push!(cur2_l,cur2)
		push!(scur1_l,scur1)
		push!(scur2_l,scur2)
		push!(q_l, qn)
		push!(p_l, pn)
		push!(mom_el, mom)
		push!(loc_cur_l,loccur)
		push!(tim,abs(i*dt))
	end

	if i % GlobalConstants.SaveData == 0
		convfile = open("$(path)_convergence_info.dat","a")		
		curfile = open("$(path)_currents.dat","a")
		popfile = open("$(path)_electrons.dat","a")
		popesfile = open("$(path)_electrons_es.dat","a")
		phfile  = open("$(path)_phonons.dat","a")
		momfile = open("$(path)_momentum_el.dat","a")
		loccurf = open("$(path)_loc_cur.dat","a")
		for j in 1:length(pop_l)
			write(convfile,"$(tim[j]) $(real(conv_l[j][1])) $(real(conv_l[j][2])) $(real(conv_l[j][3])) \n")
			write(popfile,"$(tim[j]) ")
			write(popesfile,"$(tim[j]) ")
			write(curfile,"$(tim[j]) $(cur1_l[j]) $(cur2_l[j]) $(scur1_l[j][1]) $(scur1_l[j][2]) $(scur1_l[j][3]) $(scur2_l[j][1]) $(scur2_l[j][2]) $(scur2_l[j][3]) \n")
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
		close(convfile)
		close(curfile)
		close(popfile)		
		close(popesfile)	
		close(phfile)
		close(momfile)
		close(loccurf)
		conv_l = []
		pop_l = []
		pop_es_l = []
		cur1_l = []
		cur2_l = []
		scur1_l = []
		scur2_l = []
		q_l = []
		p_l = []
		mom_el = []
		loc_cur_l = []
		tim = []

	end


	return conv_l,pop_l,pop_es_l,tim,cur1_l,cur2_l,scur1_l,scur2_l,q_l,p_l,mom_el,loc_cur_l
end

function Average(ind,path,num_avg)
	pathr = "$(path)_$(ind)"
	observables = ["currents","electrons","electrons_es","phonons","momentum_el","loc_cur"]
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
			F_avg[2:end,2:end] += file_raw[2:end,2:end] * (ind-1)
			F_avg[2:end,2:end] = F_avg[2:end,2:end] ./ ind
			rm("$(path)_$(obs)_avg.dat")
			writedlm("$(path)_$(obs)_avg.dat",F_avg) 
			rm("$(pathr)_$(obs).dat")
		end
	end
end

end
