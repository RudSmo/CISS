module Electrons
include("global_constants.jl")
using .GlobalConstants
using LinearAlgebra
export Populations,PopulationsES,DensityFEx, VNeumann, VNeumannL, MomentumElectrons, LocalCurrent, LocalSpinCoherence, LocalSpinDensity, PadeFermi, GammaMatrices1d, GammaMatrices2d, CurrentMatrix, CurrentEOM, CurrentInit_PFree, Current, DensityFExOpen, SpinCurrent

# Initial state

function DensityFEx(lx  ::Int,
                    ly  ::Int,
                    hcen::Array,
                    fil :: Int,
                    ex  :: Int)
        eigvec = eigvecs(hcen)
        teigvc = conj(transpose(eigvec))
        rho = zeros(2*lx*ly,2*lx*ly)*(1+1im*0)
        for i in 1:Int(2*lx*ly)
                for j in 1:Int(2*lx*ly)
                        rho[i,j] = sum([eigvec[i,a]*teigvc[a,j] for a in 1:(fil)])#+eigvec[i,fil+ex]*teigvc[fil+ex,j]
                end
        end
        return rho
end

# Observables

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

function LocalSpinCoherence(rho)
       sx = [0 1; 1 0.]
       sy = [0 -1im; 1im 0]
       sz = [1 0; 0 -1.]
       N = size(rho)[1]
       coh = [rho[i:i+1,i+2:i+3] + rho[i+2:i+3,i:i+1] for i in 1:2:N-3]
       return [0.5 .* [tr(sx * ci), tr(sy * ci),tr(sz * ci)] for ci in coh]
end

function LocalSpinDensity(rho)	
       	sx = [0 1; 1 0.]
       	sy = [0 -1im; 1im 0]
       	sz = [1 0; 0 -1.]
       	N = size(rho)[1]
       	dens = [rho[i:i+1,i:i+1] for i in 1:2:N]
       	return [0.5 .* [tr(sx * ci), tr(sy * ci),tr(sz * ci)] for ci in dens]	
end

# Equation of motion

function VNeumann(rho,hcen)
	return -1im*(hcen*rho-rho*hcen)
end

function VNeumannL(rho,H,PiL)
	return -1*1im*(H*rho-rho*H) + PiL + conj(transpose(PiL))
end


# Pade decomposition

function B(m)
        return 2.0*m-1
end

function PadeFermi(npade)
        diagon = fill(0.0,2*npade)
        ofdiag = [1/( sqrt( B(i)*B(i+1) ) ) for i in 1:2*npade-1]
        p = SymTridiagonal(diagon,ofdiag)
        ewp,evp = eigen(p)
        diagon = ewp
        eig_denum = diagon[npade+1:2*npade]
        eig_denum = 2.0 ./eig_denum
        ofdiag = [1/( sqrt( B(i)*B(i+1) ) ) for i in 2:2*npade-1]
        diagon = fill(0.0,2*npade-1)
        p = SymTridiagonal(diagon,ofdiag)
        ewp,evp = eigen(p)
        diagon = ewp
        ofdiag = fill(0.0,2*npade-2)
        eig_num = diagon[npade+1:2*npade-1]
        eig_num = 2.0 ./ eig_num
        kappa = [1.0 for i in 1:npade]
        for l1 in 1:npade
                for l2 in 1:npade-1
                        if l2 != l1
                                kappa[l1] = kappa[l1] * ( ( eig_num[l2]^2 - eig_denum[l1]^2 ) / ( eig_denum[l2]^2 - eig_denum[l1]^2 ) )
                        else
                                kappa[l1] = kappa[l1] * (eig_num[l2]^2 - eig_denum[l1]^2)
                        end
                end
                if l1 != npade
                        kappa[l1] = kappa[l1] / (eig_denum[npade]^2 - eig_denum[l1]^2)
                end
                kappa[l1] = kappa[l1] * 0.5 * npade * B(npade+1)
        end
        return kappa,eig_denum
end



function GammaMatrices1d(l,tL,tR)
	G1 = zeros(2*l,2*l)
	G2 = zeros(2*l,2*l)
	G1[1,1] = tL[1]
	G1[2,2] = tL[2]
	G2[2*l-1,2*l-1] = tR[1]
	G2[2*l,2*l] = tR[2]
	return G1,G2
end

function GammaMatrices2d(l,lx,tL,tR)
	G1 = zeros(2*l,2*l)
	G2 = zeros(2*l,2*l)
	G1[1,1] = tL[1]
	G1[2,2] = tL[2]
	bounds_l = [i for i in 1:2*l if i % (2*lx+1) == 0]
	bounds_r = [i for i in 1:2*l if i % (2*lx) == 0]
	for b in bounds_l
		G1[b,b] = tL[1]
		G1[b+1,b+1] = tL[2]
	end
	for b in bounds_r 
		G2[b-1,b-1] = tR[1]
		G2[b,b] = tR[2]
	end
	return G1,G2
end


function CurrentMatrix(rho,GammaL,CurrentMatrixP)
	n  = Int(sqrt(length(rho)))
	id = Matrix{Float64}(I,n,n)
	return 1/4 * (id-2*rho)*GammaL+sum(CurrentMatrixP) #* GammaL # DELETE HERE IF TYPO
end


function CurrentEOM(beta,GammaL,H,GammaT,EtaLP,ChiLP,CurrentMatrixP)
	n  = Int(sqrt(length(H)))
	id = Matrix{Float64}(I,n,n)
	return -1im*(EtaLP/beta * GammaL + (H - 1im/2 * GammaT - ChiLP * id) * CurrentMatrixP)
end

function CurrentInit_PFree(hcen,GammaT,chiLP,etaLP,beta,GammaL)
	id = Matrix{Float64}(I,length(hcen[1,:]),length(hcen[1,:]))
	inV = inv(hcen-1im*0.5*GammaT-chiLP*id)
	return etaLP/beta * inV * GammaL
end


function Current(CurrentMatrix)
	return real(tr(CurrentMatrix))
end


function DensityFExOpen(lx  ::Int,
                    ly  ::Int,
                    hcen::Array,
		    etaLP::Array,
		    chiLP::Array,
		    GammaT::Array,
		    beta ::Float64)
	id = Matrix{Float64}(I,Int(lx*ly*2),Int(lx*ly*2))
	function DensMatInv(chilp,hcen,GammaT,etalp,beta)
		inV = inv(chilp*id-hcen+0.5*1im*GammaT)
		sum = real(etalp/beta * inV)
		return sum
	end
	rho = 1/2 * id + sum([DensMatInv(chiLP[1][i],hcen,GammaT,etaLP[i],beta) for i in 1:length(chiLP[1])]) + sum([DensMatInv(chiLP[2][i],hcen,GammaT,etaLP[i],beta) for i in 1:length(chiLP[2])])
    return rho
end

function SpinCurrent(CurrentMatrix)
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
