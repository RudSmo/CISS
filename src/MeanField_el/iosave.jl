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

#function SaveOutput(cur1_c,cur2_c,pop_l,rho_i,tim,i,dt,cur1_l,cur2_l,scur1_l,scur2_l,path)
function SaveOutput(cur1_c,cur2_c,pop_l,rho_i,tim,i,dt,cur1_l,cur2_l,cur1u_l,cur1d_l,cur2u_l,cur2d_l,scur1_l,scur2_l,path)

	pop = Populations(rho_i)
	cur1 = Current(cur1_c)
	cur2 = Current(cur2_c)

	cur1u = CurrentSP(cur1_c)[1]
 	cur1d = CurrentSP(cur1_c)[2]
	cur2u = CurrentSP(cur2_c)[1]
        cur2d = CurrentSP(cur2_c)[2]

	scur1 = SpinCurrent(cur1_c)
        scur2 = SpinCurrent(cur2_c)
	if i == 1
		curfile = open("$(path)_currents.dat","w")
		popfile = open("$(path)_electrons.dat","w")
#		write(curfile,"#(1)t (2)IL (3)IR (4)sILx (5)sILy (6)sILz (7)sIRx (8)sIRy (9)sIRz\n")	
		write(curfile,"#(1)t (2)IL (3)ILu (4)ILd (5)IR (6)IRu (7)IRd (8)sILx (9)sILy (10)sILz (11)sIRx (12)sIRy (13)sIRz\n")
		write(popfile,"#(1)t ")
		k = 2
		for m in 1:length(pop)
			write(popfile,"($(k))p$(m)u ($(k+1))p$(m)d ")
			k += 2
		end
		write(popfile,"\n")

		close(curfile)
		close(popfile)
	end
	
	if i % GlobalConstants.AppendData == 0
		push!(pop_l,pop)
		push!(cur1_l,cur1)
		push!(cur2_l,cur2)
		push!(cur1u_l,cur1u)
		push!(cur1d_l,cur1d)
                push!(cur2u_l,cur2u)
		push!(cur2d_l,cur2d)
		push!(scur1_l,scur1)
                push!(scur2_l,scur2)	
		push!(tim,abs(i*dt))
	end

	if i % GlobalConstants.SaveData == 0
		curfile = open("$(path)_currents.dat","a")
		popfile = open("$(path)_electrons.dat","a")
			
		for j in 1:length(pop_l)
			write(popfile,"$(tim[j]) ")
			write(curfile,"$(tim[j]) $(cur1_l[j]) $(cur1u_l[j]) $(cur1d_l[j]) $(cur2_l[j]) $(cur2u_l[j]) $(cur2d_l[j]) $(scur1_l[j][1]) $(scur1_l[j][2]) $(scur1_l[j][3]) $(scur2_l[j][1]) $(scur2_l[j][2]) $(scur2_l[j][3]) \n")
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
		cur1u_l = []
		cur1d_l = []
                cur2u_l = []
		cur2d_l = []
		scur1_l = []
		scur2_l = []
		tim = []

	end


#	return pop_l,tim,cur1_l,cur2_l,scur1_l,scur2_l
	return pop_l,tim,cur1_l,cur2_l,cur1u_l,cur1d_l,cur2u_l,cur2d_l,scur1_l,scur2_l
end


end
