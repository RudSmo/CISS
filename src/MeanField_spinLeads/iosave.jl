module IOsave
include("electrons.jl")
include("global_constants.jl")
using .GlobalConstants
using .Electrons
using LinearAlgebra
using Dates
using DelimitedFiles
export SaveParams,SaveOutput,SaveOutputSpin,Average

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

function SaveOutput(cur1_c,cur2_c,pop_l,rho_i,qn,pn,tim,i,dt,cur1_l,cur2_l,scur1_l,scur2_l,q_l,p_l,path)

	pop = Populations(rho_i)
	cur1 = Current(cur1_c)
	cur2 = Current(cur2_c)
	scur1 = SpinCurrent(cur1_c)
	scur2 = SpinCurrent(cur2_c)

	if i == 1
		curfile = open("$(path)_currents.dat","w")
		popfile = open("$(path)_electrons.dat","w")
		phfile  = open("$(path)_phonons.dat","w")
		write(curfile,"#(1)t (2)IL (3)IR (4)sILx (5)sILy (6)sILz (7)sIRx (8)sIRy (9)sIRz\n")	
		write(popfile,"#(1)t ")
		for k in 1:length(qn)
			write(popfile,"($(k+1))p$(k)u ($(k+1))p$(k)d ")
		end
		write(popfile,"\n")
		write(phfile,"(1)t ")
		k = 2
		for i in 1:length(qn)
			write(phfile,"($(k))q$(i) ($(k+1))p$(i) ")
			k += 2
		end
		write(phfile,"\n")

		close(curfile)
		close(popfile)
		close(phfile)
	end
	
	if i % GlobalConstants.AppendData == 0
		push!(pop_l,pop)
		push!(cur1_l,cur1)
		push!(cur2_l,cur2)
		push!(scur1_l,scur1)
		push!(scur2_l,scur2)
		push!(q_l, qn)
		push!(p_l, pn)
		push!(tim,abs(i*dt))
	end

	if i % GlobalConstants.SaveData == 0
		curfile = open("$(path)_currents.dat","a")
		popfile = open("$(path)_electrons.dat","a")
		phfile  = open("$(path)_phonons.dat","a")
			
		for j in 1:length(pop_l)
			write(popfile,"$(tim[j]) ")
			write(curfile,"$(tim[j]) $(cur1_l[j]) $(cur2_l[j]) $(scur1_l[j][1]) $(scur1_l[j][2]) $(scur1_l[j][3]) $(scur2_l[j][1]) $(scur2_l[j][2]) $(scur2_l[j][3]) \n")
			write(phfile, "$(tim[j]) ")
			for k in 1:length(pop_l[j])
				write(popfile,"$(pop_l[j][k]) ")
			end
			for k in 1:length(q_l[1])
				write(phfile, "$(q_l[j][k]) $(p_l[j][k]) ")
			end
			write(popfile, "\n")
			write(phfile, "\n")
		end
		close(curfile)
		close(popfile)			
		close(phfile)
		pop_l = []
		cur1_l = []
		cur2_l = []
		scur1_l = []
		scur2_l = []
		q_l = []
		p_l = []
		tim = []

	end


	return pop_l,tim,cur1_l,cur2_l,scur1_l,scur2_l,q_l,p_l
end

#function SaveOutputSpin(cur1u_c,cur1d_c,cur2u_c,cur2d_c,pop_l,rho_i,qn,pn,tim,i,dt,cur1u_l,cur1d_l,cur2u_l,cur2d_l,scur1u_l,scur1d_l,scur2u_l,scur2d_l,q_l,p_l,path)
function SaveOutputSpin(cur1u_c,cur1d_c,cur2u_c,cur2d_c,pop_l,rho_i,qn,pn,tim,i,dt,cur1_l,cur2_l,scur1_l,scur2_l,q_l,p_l,path)

	pop = Populations(rho_i)

	cur1 = Current(cur1u_c+cur1d_c)
#	cur1u = Current(cur1u_c) 
#	cur1d = Current(cur1d_c)	

	cur2 = Current(cur2u_c+cur2d_c)
#	cur2d = Current(cur2u_c)
#	cur2d = Current(cur2d_c)
	
	scur1 = SpinCurrent(cur1u_c+cur1d_c)
	#scur1d = SpinCurrent(cur1d_c)
		
	scur2 = SpinCurrent(cur2u_c+cur2d_c)

	#scur2d = SpinCurrent(cur2d_c)
		

	if i == 1
		curfile = open("$(path)_currents.dat","w")
		popfile = open("$(path)_electrons.dat","w")
		phfile  = open("$(path)_phonons.dat","w")
		#write(curfile,"#(1)t (2)ILu (3)ILd (4)IRu (5)IRd (6)sILux (7)sILuy (8)sILuz (9)sILdx (10)sILdy (11)sILdz (12)sIRux (13)sIRuy (14)sIRuz (15)sIRdx (16)sIRdy (17)sIRdz\n")	
		write(curfile,"#(1)t (2)IL (3)IR (4)sILx (5)sILy (6)sILz (7)sIRx (8)sIRy (9)sIRz\n")	
		write(popfile,"#(1)t ")
		for k in 1:length(qn)
			write(popfile,"($(k+1))p$(k)u ($(k+1))p$(k)d ")
		end
		write(popfile,"\n")
		write(phfile,"(1)t ")
		k = 2
		for i in 1:length(qn)
			write(phfile,"($(k))q$(i) ($(k+1))p$(i) ")
			k += 2
		end
		write(phfile,"\n")

		close(curfile)
		close(popfile)
		close(phfile)
	end
	
	if i % GlobalConstants.AppendData == 0
		push!(pop_l,pop)
		push!(cur1_l,cur1)
		#push!(cur1u_l,cur1u)
		#push!(cur1d_l,cur1d)
		
		push!(cur2_l,cur2)
	#	push!(cur2u_l,cur2u)
	#	push!(cur2d_l,cur2d)
		
		push!(scur1_l,scur1)	
#		push!(scur1u_l,scur1u)
#		push!(scur1d_l,scur1d)
	
		push!(scur2_l,scur2)
	#	push!(scur2u_l,scur2u)
	#	push!(scur2d_l,scur2d)
			
		push!(q_l, qn)
		push!(p_l, pn)
		push!(tim,abs(i*dt))
	end

	if i % GlobalConstants.SaveData == 0
		curfile = open("$(path)_currents.dat","a")
		popfile = open("$(path)_electrons.dat","a")
		phfile  = open("$(path)_phonons.dat","a")
			
		for j in 1:length(pop_l)
			write(popfile,"$(tim[j]) ")
			#write(curfile,"$(tim[j]) $(cur1u_l[j]) $(cur1d_l[j]) $(cur2u_l[j]) $(cur2d_l[j]) $(scur1u_l[j][1]) $(scur1u_l[j][2]) $(scur1u_l[j][3]) $(scur1d_l[j][1]) $(scur1d_l[j][2]) $(scur1d_l[j][3]) $(scur2u_l[j][1]) $(scur2u_l[j][2]) $(scur2u_l[j][3]) $(scur2d_l[j][1]) $(scur2d_l[j][2]) $(scur2d_l[j][3])\n")
			write(curfile,"$(tim[j]) $(cur1_l[j]) $(cur2_l[j]) $(scur1_l[j][1]) $(scur1_l[j][2]) $(scur1_l[j][3]) $(scur2_l[j][1]) $(scur2_l[j][2]) $(scur2_l[j][3]) \n")
				
			write(phfile, "$(tim[j]) ")
			for k in 1:length(pop_l[j])
				write(popfile,"$(pop_l[j][k]) ")
			end
			for k in 1:length(q_l[1])
				write(phfile, "$(q_l[j][k]) $(p_l[j][k]) ")
			end
			write(popfile, "\n")
			write(phfile, "\n")
		end
		close(curfile)
		close(popfile)			
		close(phfile)

		pop_l = []
		
		cur1_l = []
		#cur1u_l = []
		#cur1d_l = []
	
		cur2_l = []
	#	cur2u_l = []
	#	cur2d_l = []
		
		scur1_l = []
	#	scur1u_l = []
	#	scur1d_l = []
	#		
		scur2_l = []
	#	scur2u_l = []
	#	scur2d_l = []
			
		q_l = []
		p_l = []
		tim = []

	end


	#return pop_l,tim,cur1u_l,cur1d_l,cur2u_l,cur2d_l,scur1u_l,scur1d_l,scur2u_l,scur2d_l,q_l,p_l
	return pop_l,tim,cur1_l,cur2_l,scur1_l,scur2_l,q_l,p_l

end


function Average(ind,path,num_avg)
	pathr = "$(path)_$(ind)"
	observables = ["currents","electrons","phonons"]
	rm("$(pathr)_final_Quantities")
	if ind == 1
		for obs in observables
			file_avg = readdlm("$(pathr)_$(obs).dat")
			file_avg[2:end,2:end] = file_avg[2:end,2:end] ./ num_avg
			writedlm("$(path)_$(obs)_avg.dat",file_avg) 
			rm("$(pathr)_$(obs).dat")
		end
	else
		for obs in observables
			file_avg = readdlm("$(path)_$(obs)_avg.dat")
			file_raw = readdlm("$(pathr)_$(obs).dat")

			F_avg    = copy(file_avg)
			F_avg[2:end,2:end] += file_raw[2:end,2:end] ./ num_avg
			rm("$(path)_$(obs)_avg.dat")
			writedlm("$(path)_$(obs)_avg.dat",F_avg) 
			rm("$(pathr)_$(obs).dat")
		end
	end

end

end
