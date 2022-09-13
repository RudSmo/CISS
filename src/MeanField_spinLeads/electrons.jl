module Electrons
include("global_constants.jl")
using .GlobalConstants
using LinearAlgebra
export Populations,DensityFEx, VNeumann, VNeumannL, PadeFermi, GammaMatrices1d, GammaMatrices1dSpin, GammaMatrices2d, CurrentMatrix, CurrentEOM, CurrentInit_PFree, Current, DensityFExOpen, SpinCurrent

function Populations(rho)
	return [real(rho[i,i]) for i in 1:size(rho)[1]]
end

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


function VNeumann(rho,hcen)
	return -1im*(hcen*rho-rho*hcen)
end

function VNeumannL(rho,H,PiL)
	return -1*1im*(H*rho-rho*H) + PiL + conj(transpose(PiL))
end

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
function GammaMatrices1dSpin(l,tL,tR)
	G1u = zeros(2*l,2*l)
	G1d = zeros(2*l,2*l)
	G2u = zeros(2*l,2*l)
	G2d = zeros(2*l,2*l)
		
	G1u[1,1] = tL[1,1]
	G1u[1,2] = tL[1,2]
	G1d[2,1] = tL[2,1]
	G1d[2,2] = tL[2,2]
	
	G2u[2*l-1,2*l-1] = tR[1,1]
	G2u[2*l-1,2*l] 	 = tR[1,2]
	G2d[2*l,2*l-1] 	 = tR[2,1]
	G2d[2*l,2*l] 	 = tR[2,2]


	return G1u,G1d,G2u,G2d
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
                    ly  	::Int,
                    hcen	::Array,
		    		etaLP	::Array,
		    		chiLP	::Array,
		    		GammaT	::Array,
		    		beta 	::Array)
	id = Matrix{Float64}(I,Int(lx*ly*2),Int(lx*ly*2))
	function DensMatInv(chilp,hcen,GammaT,etalp,beta)
		inV = inv(chilp*id-hcen+0.5*1im*GammaT)
		sum = real(etalp/beta * inV)
		return sum
	end
	rho = 1/2 * id + sum([sum([DensMatInv(chiLP[j][i],hcen,GammaT[j],etaLP[i],beta[j]) for i in 1:length(chiLP[j])]) for j in 1:length(chiLP)]) 
	#+ 
#	                 sum([DensMatInv(chiLP[2][i],hcen,GammaT,etaLP[i],beta) for i in 1:length(chiLP[2])]) +
#					 sum([DensMatInv(chiLP[3][i],hcen,GammaT,etaLP[i],beta) for i in 1:length(chiLP[3])]) + 
#	                 sum([DensMatInv(chiLP[4][i],hcen,GammaT,etaLP[i],beta) for i in 1:length(chiLP[4])])
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
