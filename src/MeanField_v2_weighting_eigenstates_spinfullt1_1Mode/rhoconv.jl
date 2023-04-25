module RhoConvert

include("global_constants.jl")
using .GlobalConstants
using LinearAlgebra
export RhoConv,RhoConvMulti,RhoConvTotal

function RhoConv(rho,pos)
	""" Takes a density matrix, partitions it into 2x2 spin blocks mutliplied by the corresponding 
	 Pauli matrices (sintx,sinty,sintz). Takes the trace of those, multiplied with a factor of 0.5.
	 Here, the spin polarization of the conduction electron system is computed.""" 
	 #S = []
	 #for j in 1:length(pos)
	 i      = 2*pos-1
	 rhoi   = [rho[i,i] rho[i,i+1];rho[i+1,i] rho[i+1,i+1]]
	 sp     = [0.5*tr(rhoi*GlobalConstants.s_x),0.5*tr(rhoi*GlobalConstants.s_y),0.5*tr(rhoi*GlobalConstants.s_z)]
	 return sp
end

function RhoConvMulti(rho,pos)
	""" Takes a density matrix, partitions it into 2x2 spin blocks mutliplied by the corresponding 
	 Pauli matrices (sintx,sinty,sintz). Takes the trace of those, multiplied with a factor of 0.5.
	 Here, the spin polarization of the conduction electron system is computed.""" 
	 #S = []
	 #for j in 1:length(pos)
	 S = []
	 for j in pos
		 i      = 2*j-1
		 rhoi   = [rho[i,i] rho[i,i+1];rho[i+1,i] rho[i+1,i+1]]
	 	 sp     = [0.5*tr(rhoi*GlobalConstants.s_x),0.5*tr(rhoi*GlobalConstants.s_y),0.5*tr(rhoi*GlobalConstants.s_z)]
		 push!(S,sp)
	 end
	 return S
end

function RhoConvTotal(rho)
	rhoi = [[rho[i,i] rho[i,i+1];rho[i+1,i] rho[i+1,i+1]] for i in 1:2:Int(sqrt(length(rho)))]
	sp = [[0.5*tr(rhoi[i]*GlobalConstants.s_x),0.5*tr(rhoi[i]*GlobalConstants.s_y),0.5*tr(rhoi[i]*GlobalConstants.s_z)] for i in 1:length(rhoi)]
	return sp
end


end

