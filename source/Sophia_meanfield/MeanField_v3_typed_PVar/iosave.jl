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
	#loc_cur_l::Array{Float64},
	q_l::Array{Float64},
	p_l::Array{Float64},
	path::String,
	vind::Int64)::Tuple{Array{Float64},Array{Float64},Array{Float64},Array{Float64},Array{Float64},Array{Float64},Array{Float64},Array{Float64},Array{Float64},Array{Float64},Int64}#Array{Float64}

if i == 1
convfile= open("$(path)_convergence_info.dat","w")
curfile = open("$(path)_currents.dat","w")
popfile = open("$(path)_electrons.dat","w")
popesfile = open("$(path)_electrons_es.dat","w")
phfile  = open("$(path)_phonons.dat","w")

#loccurf = open("$(path)_loc_cur.dat","w")

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



#=write(loccurf,"#(1)t")
for i in 1:HamPar.N-1
	write(loccurf,"($(i+1))cur($(i),$(i+1)) ")
end
write(loccurf,"\n")=#

close(convfile)
close(curfile)
close(popfile)
close(popesfile)
close(phfile)
#close(loccurf)
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
#loc_cur_l[:,vind] = LocalCurrent(rho_i,hcen_i)
tim[vind] = abs(i*dt)
vind += 1
end

if i % SaveData == 0
convfile = open("$(path)_convergence_info.dat","a")		
curfile = open("$(path)_currents.dat","a")
popfile = open("$(path)_electrons.dat","a")
popesfile = open("$(path)_electrons_es.dat","a")
phfile  = open("$(path)_phonons.dat","a")

#loccurf = open("$(path)_loc_cur.dat","a")
for j in 1:length(tim)
	write(convfile,"$(tim[j]) $(real(conv_l[1,j])) $(real(conv_l[2,j])) $(real(conv_l[3,j])) \n")
	write(popfile,"$(tim[j]) ")
	write(popesfile,"$(tim[j]) ")
	write(curfile,"$(tim[j]) $(cur1_l[j]) $(cur2_l[j]) $(scur1_l[1,j]) $(scur1_l[2,j]) $(scur1_l[3,j]) $(scur2_l[1,j]) $(scur2_l[2,j]) $(scur2_l[3,j]) \n")
	write(phfile, "$(tim[j]) ")
	
	#write(loccurf, "$(tim[j]) ")

		
	for k in 1:length(pop_l[:,j])
		write(popfile,"$(pop_l[k,j]) ")
		write(popesfile,"$(pop_es_l[k,j]) ")
	end
	for k in 1:length(q_l[:,j])
		write(phfile, "$(q_l[k,j]) $(p_l[k,j]) ")
	end
	#=
	for k in 1:length(loc_cur_l[:,j])
		write(loccurf, "$(loc_cur_l[k,j]) ")
	end=#

	write(popfile, "\n")
	write(popesfile, "\n")
	write(phfile, "\n")
	
	#write(loccurf, "\n")
end
close(convfile)
close(curfile)
close(popfile)		
close(popesfile)	
close(phfile)
#close(loccurf)
vind = 1
end

return conv_l,pop_l,pop_es_l,tim,cur1_l,cur2_l,scur1_l,scur2_l,q_l,p_l,vind#,loc_cur_l
end
function SaveOutput_h5(HamPar::Any,
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
		    path::String,
		    vind::Int64,
			)::Tuple{Array{Float64},Array{Float64},Array{Float64},Array{Float64},Array{Float64},Array{Float64},Array{Float64},Array{Float64},Array{Float64},Array{Float64},Array{Float64},Int64}

	if i == 1
		
		fil = h5open("$(path)_final_Quantities.h5","r+")
		fil["convfile"] = ["#(1)t (2)trace (3)sys_energy (4)dif\n"] # write schreibt einträge ins doc. Hier die Kopfzeie mit Benneung der Spalten
		fil["curfile"] = ["#(1)t (2)IL (3)IR (4)sILx (5)sILy (6)sILz (7)sIRx (8)sIRy (9)sIRz\n"]	# der \n" befehl ist ein newline
		fil["popfile"]=["#(1)t "]
		fil["popesfile"] =["#(1)t "]
		fil["phfile"] =[]
	

		#what actually does this???????????????????????????
		for k in 1:HamPar.N
			write(fil["popfile"],"($(k+1))p$(k)u ($(k+1))p$(k)d ")
		end
		for k in 1:Int(2*HamPar.N)
                        write(fil["popesfile"],"($(k+1))S$(k) ")
                end	
		
		write(fil["popfile"],"\n")
		write(fil["popesfile"],"\n")
		write(fil["phfile"],"#(1)t ")
		k = 2
		for i in 1:HamPar.N
			write(fil["phfile"],"($(k))q$(i) ($(k+1))p$(i) ")
			k += 2
		end
		write(fil["phfile"],"\n")

		
		k = 2
		for i in 1:HamPar.N-1
			
			k += 2
		end	
		
		close(fil)
	end
	
	if i % AppendData == 0 # zwischenspeichert Ergebnisse für  30 Zeitschritte
# several temporary lists
		fil = h5open("$(path)_final_Quantities.h5","r+")
		fil["conv_l"][:,vind] = ConvergenceInfo(HamPar,rho_i,hcen_i)
		fil["pop_l"][:,vind] = Populations(rho_i)
		fil["pop_es_l"][:,vind]= PopulationsES(rho_i,hcen_i)
		fil["pop_l"][:, vind] = Populations(rho_i)
		fil["cur1_l"][:, vind] = Current(cur1_c)
		fil["cur2_l"][:, vind] = Current(cur2_c)
		fil["scur1_l"][:, :, vind] = SpinCurrent(cur1_c)
		fil["scur2_l"][:, :, vind] = SpinCurrent(cur2_c)
		fil["q_l"][:, :, vind] = HamPar.X
		fil["p_l"][:, :, vind] = HamPar.P
		#fil["loc_cur_l"][:, :, vind] = LocalCurrent(rho_i, hcen_i)
		fil["tim"][vind] = abs(i*dt)
		vind += 1 
		close(fil)
	end

	if i % SaveData == 0 #saved die Zwischenspeicherung der letzten 30 Zeitschritte
		fil = h5open("$(path)_final_Quantities.h5","r+")
		
		
		for j in 1:length(tim)
			write(fil["convfile"],"$(fil["tim"][j]) $(real(fil["conv_l"][1,j])) $(real(fil["conv_l"][2,j])) $(real(fil["conv_l"][3,j])) \n")
			write(fil["popfile"],"$(tim[j]) ") 
			write(fil["popesfile"],"$(tim[j]) ")
			write(fil["curfile"],"$(tim[j]) $(fil["cur1_l"][j]) $(fil["cur2_l"][j]) $(fil["scur1_l"][1, j]) $(fil["scur1_l"][2, j]) $(fil["scur1_l"][3, j]) $(fil["scur2_l"][1, j]) $(fil["scur2_l"][2, j]) $(fil["scur2_l"][3, j]) \n")
			write(fil["phfile"], "$(tim[j]) ")
			#write(fil["loccurf"], "$(tim[j]) ")

			for k in 1:length(fil["pop_l"][:, j])
				write(fil["popfile"],"$(fil["pop_l"][k, j]) ")  
				write(fil["popesfile"],"$(fil["pop_es_l"][k, j]) ")
			end
			for k in 1:length(fil["q_l"][:, j])
				write(fil["phfile"], "$(fil["q_l"][k, j]) $(fil["p_l"][k, j]) ")
			end

			#for k in 1:length(fil["loc_cur_l"][:, j])
			#	write(fil["loccurf"], "$(fil["loc_cur_l"][k, j]) ")
			#end

			write(fil["popfile"], "\n")
			write(fil["popesfile"], "\n")
			write(fil["phfile"], "\n")
			#write(fil["loccurf"], "\n")
		end

		close(fil)
		vind = 1  
	end

	return fil["conv_l"], fil["pop_l"], fil["pop_es_l"], fil["pop_l"], fil["cur1_l"],fil["cur2_l"], fil["scur1_l"],fil["scur2_l"], fil["q_l"],	fil["p_l"], fil["tim"]
	end

function Average(ind,path,num_avg) # mittelt die Werte von verschiedenen Runs. Die Startwerte werden nach dem Zufall initialisiert, daher ist das average wichtig für ein statistisch genaues Ergebnis
	pathr = "$(path)_$(ind)"
	observables = ["currents","electrons","electrons_es","phonons","loc_cur","convergence_info"]
	rm("$(pathr)_final_Quantities")
	if ind == 1
		for obs in observables
			file_avg = readdlm("$(pathr)_$(obs).dat")
			writedlm("$(path)_$(obs)_avg.dat",file_avg) 
			rm("$(pathr)_$(obs).dat")
		end
	elseif ind < num_avg
	
		for obs in observables
			file_avg = readdlm("$(path)_$(obs)_avg.dat")
			file_raw = readdlm("$(pathr)_$(obs).dat")

			F_avg    = copy(file_avg)
			for k in 2:length(F_avg[2,:])
				F_avg[2:end,k] .+= file_raw[2:end,k]
			end
			rm("$(path)_$(obs)_avg.dat")
			writedlm("$(path)_$(obs)_avg.dat",F_avg) 
			rm("$(pathr)_$(obs).dat")
		end
	else
		for obs in observables
                        file_avg = readdlm("$(path)_$(obs)_avg.dat")
                        file_raw = readdlm("$(pathr)_$(obs).dat")

                        F_avg    = copy(file_avg)
			for k in 2:length(F_avg[2,:])
				F_avg[2:end,k] .+= file_raw[2:end,k]
				F_avg[2:end,k] ./= num_avg
			end
                        rm("$(path)_$(obs)_avg.dat")
                        writedlm("$(path)_$(obs)_avg.dat",F_avg)
                        rm("$(pathr)_$(obs).dat")
                end
	end
end

