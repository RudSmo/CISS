module SaveOutput

include("iosave.jl")
include("electrons.jl")
include("global_constants.jl")
using .IOsave
using .Electrons
using .GlobalConstants

export SaveMultiClosed,SaveMultiOpen


function SaveMultiClosed(rho_l,rho_i,scd_l,scd_i,mag_l,mag_i,tim,i,dt,writeMode,pos,path,L)
	if i == 1
		push!(rho_l,rho_i)
		push!(scd_l,scd_i)
		push!(mag_l,real.(mag_i))
		push!(tim,i*dt)
		if writeMode == "minimal"
			IOsave.MultiIntegratedMinimal(mag_l,scd_l,pos,path,tim,"w")
		end
		if writeMode == "populations"
			IOsave.MultiIntegratedPopulations(mag_l,scd_l,rho_l,pos,L,path,tim,"w")
		end
		rho_l = []
		scd_l = []
		mag_l = []
		tim   = []
	
	end
	if i%GlobalConstants.AppendData == 0
		
		push!(rho_l,rho_i)
		push!(scd_l,scd_i)
		push!(mag_l,mag_i)
		push!(tim,i*dt)
	end
	if i%GlobalConstants.SaveData == 0 && i != 1
		if writeMode == "minimal"
			IOsave.MultiIntegratedMinimal(mag_l,scd_l,pos,path,tim,"a")
		end
		if writeMode == "populations"
			IOsave.MultiIntegratedPopulations(mag_l,scd_l,rho_l,pos,L,path,tim,"a")
		end


		rho_l = []
		scd_l = []
		mag_l = []
		tim   = []
		
	end
	return rho_l,scd_l,mag_l,tim
end

function SaveMultiOpen(cur1_c,cur2_c,rho_l,rho_i,scd_l,scd_i,mag_l,mag_i,pos,tim,i,dt,cur1_l,cur2_l,scur1_l,scur2_l,L,path,writeMode)
	if i == 1
		current_1= Current(cur1_c)
		current_2= Current(cur2_c)
		scurrent_1= SpinCurrent(cur1_c)
		scurrent_2= SpinCurrent(cur2_c)
				
		push!(rho_l,rho_i)
		push!(scd_l,scd_i)
		push!(mag_l,real.(mag_i))
		push!(tim,i*dt)
		push!(cur1_l,current_1)
		push!(cur2_l,current_2)
		push!(scur1_l,scurrent_1)
		push!(scur2_l,scurrent_2)
		if writeMode == "minimal"
			IOsave.MultiOpenMinimal(mag_l,scd_l,pos,cur1_l,cur2_l,scur1_l,scur2_l,path,tim,"w")
		end
		if writeMode == "populations"
			IOsave.MultiOpenPopulations(mag_l,scd_l,rho_l,pos,L,cur1_l,cur2_l,scur1_l,scur2_l,path,tim,"w")
		end

		rho_l = []
		scd_l = []
		mag_l = []
		tim   = []
		cur1_l = []
		cur2_l = []
		scur1_l= []
		scur2_l= []
	end
	if i%GlobalConstants.AppendData == 0
		current_1= Current(cur1_c)
		current_2= Current(cur2_c)
		scurrent_1= SpinCurrent(cur1_c)
		scurrent_2= SpinCurrent(cur2_c)

		push!(rho_l,rho_i)
		push!(scd_l,scd_i)
		push!(mag_l,mag_i)
		push!(tim,i*dt)
		push!(cur1_l,current_1)
		push!(cur2_l,current_2)
		push!(scur1_l,scurrent_1)
		push!(scur2_l,scurrent_2)
	end
	if i%GlobalConstants.SaveData == 0 && i != 1
		if writeMode == "minimal"
			IOsave.MultiOpenMinimal(mag_l,scd_l,pos,cur1_l,cur2_l,scur1_l,scur2_l,path,tim,"a")
		end
		if writeMode == "populations"
			IOsave.MultiOpenPopulations(mag_l,scd_l,rho_l,pos,L,cur1_l,cur2_l,scur1_l,scur2_l,path,tim,"a")
		end
	
		rho_l = []
		scd_l = []
		mag_l = []
		tim   = []
		cur1_l = []
		cur2_l = []
		scur1_l= []
		scur2_l= []

	end
	return rho_l,scd_l,mag_l,tim,cur1_l,cur2_l,scur1_l,scur2_l
end

end
