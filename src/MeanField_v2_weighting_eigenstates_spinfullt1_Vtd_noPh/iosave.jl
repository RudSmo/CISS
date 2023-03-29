module IOsave
include("electrons.jl")
include("global_constants.jl")
using .GlobalConstants
using .Electrons
using LinearAlgebra
using Dates
using DelimitedFiles
export SaveParams,SaveOutput

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

function SaveOutput(cur1_c,cur2_c,pop_l,pop_es_l,rho_i,hcen_i,tim,i,dt,cur1_l,cur2_l,scur1_l,scur2_l,loc_cur_l,mom_el,path)
	
	if i == 1
		nSys = size(hcen_i)[1]
		N = Int(nSys/2)
		curfile = open("$(path)_currents.dat","w")
		popfile = open("$(path)_electrons.dat","w")
		popesfile = open("$(path)_electrons_es.dat","w")
		momentumf= open("$(path)_momentum_el.dat","w")
		loccurf = open("$(path)_loc_cur.dat","w")

		write(curfile,"#(1)t (2)IL (3)IR (4)sILx (5)sILy (6)sILz (7)sIRx (8)sIRy (9)sIRz\n")	
		write(popfile,"#(1)t ")
		write(popesfile,"#(1)t ")
		for k in 1:N
			write(popfile,"($(k+1))p$(k)u ($(k+1))p$(k)d ")
		end
		for k in 1:nSys
                        write(popesfile,"($(k+1))S$(k) ")
                end	
		
		write(popfile,"\n")
		write(popesfile,"\n")

		write(momentumf,"#(1)t ")
		k = 2
		for i in 1:N-1
			write(momentumf,"($(k))p$(i)u ($(k+1))p$(i)d ")
			k += 2
		end	
		write(momentumf,"\n")

		write(loccurf,"#(1)t")
		for i in 1:N-1
			write(loccurf,"($(i+1))cur($(i),$(i+1)) ")
		end
		write(loccurf,"\n")

		close(curfile)
		close(popfile)
		close(popesfile)
		close(momentumf)
		close(loccurf)
	end
	
	if i % GlobalConstants.AppendData == 0

		pop = Populations(rho_i)
		pop_es = PopulationsES(rho_i,hcen_i)
		mom = MomentumElectrons(rho_i)
		cur1 = Current(cur1_c)
		cur2 = Current(cur2_c)
		scur1 = SpinCurrent(cur1_c)
		scur2 = SpinCurrent(cur2_c)
		loccur = LocalCurrent(rho_i,hcen_i)
		push!(pop_l,pop)
		push!(pop_es_l,pop_es)
		push!(cur1_l,cur1)
		push!(cur2_l,cur2)
		push!(scur1_l,scur1)
		push!(scur2_l,scur2)
		push!(mom_el, mom)
		push!(loc_cur_l,loccur)
		push!(tim,abs(i*dt))
	end

	if i % GlobalConstants.SaveData == 0
		curfile = open("$(path)_currents.dat","a")
		popfile = open("$(path)_electrons.dat","a")
		popesfile = open("$(path)_electrons_es.dat","a")
		momfile = open("$(path)_momentum_el.dat","a")
		loccurf = open("$(path)_loc_cur.dat","a")
		for j in 1:length(pop_l)
			write(popfile,"$(tim[j]) ")
			write(popesfile,"$(tim[j]) ")
			write(curfile,"$(tim[j]) $(cur1_l[j]) $(cur2_l[j]) $(scur1_l[j][1]) $(scur1_l[j][2]) $(scur1_l[j][3]) $(scur2_l[j][1]) $(scur2_l[j][2]) $(scur2_l[j][3]) \n")
			write(momfile, "$(tim[j]) ")
			write(loccurf, "$(tim[j]) ")

				
			for k in 1:length(pop_l[j])
				write(popfile,"$(pop_l[j][k]) ")
				write(popesfile,"$(pop_es_l[j][k]) ")
			end
			for k in 1:length(mom_el[j])
				write(momfile, "$(mom_el[j][k]) ")
			end
			for k in 1:length(loc_cur_l[j])
				write(loccurf, "$(loc_cur_l[j][k]) ")
			end

			write(popfile, "\n")
			write(popesfile, "\n")
			write(momfile, "\n")
			write(loccurf, "\n")
		end
		close(curfile)
		close(popfile)		
		close(popesfile)	
		close(momfile)
		close(loccurf)
		pop_l = []
		pop_es_l = []
		cur1_l = []
		cur2_l = []
		scur1_l = []
		scur2_l = []
		mom_el = []
		loc_cur_l = []
		tim = []

	end


	return pop_l,pop_es_l,tim,cur1_l,cur2_l,scur1_l,scur2_l,mom_el,loc_cur_l
end


end
