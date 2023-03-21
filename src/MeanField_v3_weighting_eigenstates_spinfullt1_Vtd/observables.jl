module Observables
include("global_constants.jl")
using .GlobalConstants
using LinearAlgebra

export Populations, PopulationsES, MomentumElectrons, LocalCurrent, Current, SpinCurrent

function Populations(rho)
	return [real(rho[i,i]) for i in 1:size(rho)[1]]
end

function PopulationsES(rho,H)
	ev = eigvecs(H)
	rho_t = ev' * rho * ev
        return [real(rho_t[i,i]) for i in 1:size(rho_t)[1]]
end

function MomentumElectrons(rho)
	Q = zeros(size(rho))*(1+0*1im)
	N = size(rho)[1]
	for i in 1:N-2
		for j in 1:N
			Q[i,j] = -1im * (rho[i+2,j]-rho[i,j])
		end
	end
	return [real(abs(Q[k,k])) for k in 1:N-2] 
end

function LocalCurrent(rho,H)
	N = Int(size(rho)[1]/2)
	loccur = []
	for i in 1:N-1
		k1 = 2*i-1
		k2 = 2*(i+1)-1
		push!(loccur, 1im/2 * tr(H[k1:k1+1,k2:k2+1]*rho[k2:k2+1,k1:k1+1]-rho[k1:k1+1,k2:k2+1]*H[k2:k2+1,k1:k1+1]))		
	end
	return real.(loccur)
end

function Current(CurrentMatrix)
	# Consider changed dimension of CurrentMatrix
	return 2*real(tr(CurrentMatrix))
end

function SpinCurrent(CurrentMatrix)
	# Consider changed dimension of CurrentMatrix
	n = length(CurrentMatrix[:,1])
	a_sigmax = zeros(n,n)*(1+0*1im)
	a_sigmay = zeros(n,n)*(1+0*1im)
	a_sigmaz = zeros(n,n)*(1+0*1im)
	for i in 1:2:n
		a_sigmax[i,i+1] = GlobalConstants.s_x[1,2]
		a_sigmax[i+1,i] = GlobalConstants.s_x[2,1]		
		a_sigmay[i,i+1] = GlobalConstants.s_y[1,2]
		a_sigmay[i+1,i] = GlobalConstants.s_y[2,1]
		a_sigmaz[i,i] 	= GlobalConstants.s_z[1,1]
		a_sigmaz[i+1,i+1] = GlobalConstants.s_z[2,2]
	end
	return real.([tr(CurrentMatrix*a_sigmax),tr(CurrentMatrix*a_sigmay),tr(CurrentMatrix*a_sigmaz)])
end

end
