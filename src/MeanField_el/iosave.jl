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

function SaveOutput(cur1_c,cur2_c,pop_l,rho_i,tim,i,dt,cur1_l,cur2_l,path)

	pop = Populations(rho_i)
	cur1 = Current(cur1_c)
	cur2 = Current(cur2_c)

	if i == 1
		curfile = open("$(path)_currents.dat","w")
		popfile = open("$(path)_electrons.dat","w")
		write(curfile,"#(1)t (2)IL (3)IR \n")	
		write(popfile,"#(1)t ")
		for k in 1:length(pop)
			write(popfile,"($(k+1))p$(k) ")
		end
		write(popfile,"\n")

		close(curfile)
		close(popfile)
	end
	
	if i % GlobalConstants.AppendData == 0
		push!(pop_l,pop)
		push!(cur1_l,cur1)
		push!(cur2_l,cur2)
		push!(tim,abs(i*dt))
	end

	if i % GlobalConstants.SaveData == 0
		curfile = open("$(path)_currents.dat","a")
		popfile = open("$(path)_electrons.dat","a")
			
		for j in 1:length(pop_l)
			write(popfile,"$(tim[j]) ")
			write(curfile,"$(tim[j]) $(cur1_l[j]) $(cur2_l[j])  \n")
			for k in 1:length(pop_l[j])
				write(popfile,"$(pop_l[j][k]) ")
			end
			write(popfile, "\n")
		end
		close(curfile)
		close(popfile)			
		pop_l = []
		cur1_l = []
		cur2_l = []
		tim = []

	end


	return pop_l,tim,cur1_l,cur2_l
end


end
