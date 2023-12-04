include("global_constants.jl")
include("hamiltonian.jl")
include("classical_phonons.jl")
using LinearAlgebra
# Initial state

function DensityFEx(lx  ::Int,
                    ly  ::Int,
                    hcen::Array,
                    fil :: Int,
                    ex  :: Int)::Array
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
function ConvergenceInfo(HamPar::Any,# prüft Energieerhaltung
			rho::Array,
			H::Array)::Array{Float64}
	system_energy::ComplexF64 = tr(rho*H)
	dwork = zeros(size(H)) * (1+0*1im)
	denergy_vib::ComplexF64 = 0.0+0*1im
	for i in 1:HamPar.N
		system_energy += 0.5*(HamPar.P[i]^2 + HamPar.X[i]^2 + 0.5)
		HF_m = FranssonHamiltonian_Mode_v2(HamPar,i)
		dwork += HamPar.P[i]*HF_m
		derv = tr(@views rho*HF_m)
		denergy_vib += HamPar.X[i]*X_RHS(HamPar.P[i]) + HamPar.P[i]*P_RHS(derv,HamPar.X[i])
	end
	rd = rho * dwork
	dwork_elvib = 0.0 + 1im*0.0
	trace = 0.0 + 1im*0.0
	for j in 1:size(rho)[1]
		dwork_elvib += rd[j,j]
		trace += rho[j,j]
	end
	dif = dwork_elvib + denergy_vib
	return [real(trace),real(system_energy),real(dif)]
end


function Populations(rho::Array)::Array
	return [real(rho[i,i]) for i in 1:size(rho)[1]] # diagonaleintrage realteil
end

function PopulationsES(rho::Array,H::Array)::Array # populatiobns im energieeigenraum
	ev = eigvecs(H)
	@views rho_t = ev' * rho * ev # transformiere dichtematrix mit unitärer Energieeigenraum matrix
        return [real(rho_t[i,i]) for i in 1:size(rho_t)[1]]
end
# mom_el löschen

function LocalCurrent(rho::Array,H::Array)::Array
	N = Int(size(rho)[1]/2)
	loccur = []
	for i in 1:N-1
		k1 = 2*i-1
		k2 = 2*(i+1)-1
		push!(loccur, 1im/2 * tr(H[k1:k1+1,k2:k2+1]*rho[k2:k2+1,k1:k1+1] .- rho[k1:k1+1,k2:k2+1]*H[k2:k2+1,k1:k1+1]))		
	end
	return real.(loccur)
end

# Equation of motion

function VNeumann(rho::Array,hcen::Array)::Array
	return -1im*(hcen*rho-rho*hcen)
end

function VNeumannL(rho::Array,H::Array,PiL::Array)::Array
	return -1*1im*(H*rho-rho*H) + PiL + conj(transpose(PiL))
end

function VNeumannL_v2(rho::Array,H::Array,B::Array)::Array
	return -1im .* ( H*rho .- rho*H) .+ B
end

# Pade decomposition of  Fermi function

function B(m::Int64)::Float64
        return 2.0*m-1
end

function PadeFermi(npade::Int64)::Tuple
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



function GammaMatrices1d(l::Int64,tL::Array,tR::Array)::Tuple
	G1 = zeros(2*l,2*l)
	G2 = zeros(2*l,2*l)
	G1[1,1] = tL[1]
	G1[2,2] = tL[2]
	G2[2*l-1,2*l-1] = tR[1]
	G2[2*l,2*l] = tR[2]
	return G1,G2
end

function GammaMatrices2d(l::Int64,lx::Int64,tL::Array,tR::Array)::Tuple
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


function CurrentMatrix(rho::Array,GammaL::Array,CurrentMatrixP::Array)::Array
	n  = Int(sqrt(length(rho)))
	id = Matrix{Float64}(I,n,n)
	return 1/4 * (id-2*rho)*GammaL+sum(CurrentMatrixP,dims=3)[:,:,1]
end


function CurrentEOM(beta::Float64,GammaL::Array,H::Array,GammaT::Array,EtaLP::Float64,ChiLP::ComplexF64,CurrentMatrixP::Array)::Array
	n  = Int(sqrt(length(H)))
	id = Matrix{Float64}(I,n,n)
	return -1im*(EtaLP/beta * GammaL + (H - 1im/2 * GammaT - ChiLP * id) * CurrentMatrixP)
end

function CurrentEOM_v2(beta::Float64,GammaL::Array,H::Array,GammaT::Array,EtaLP::Float64,ChiLP::ComplexF64)::Tuple
        n  = Int(sqrt(length(H)))
        id = Matrix{Float64}(I,n,n)
        return EtaLP/beta .* GammaL, (H .- 1im/2 .* GammaT .- ChiLP .* id)
end

function CurrentEOM_RK(A::Array,B::Array,CurrentMatrixP::Array)::Array
        return -1im .* (A .+ B*CurrentMatrixP)
end

function CurrentInit_PFree(hcen::Array,GammaT::Array,chiLP::ComplexF64,etaLP::Float64,beta::Float64,GammaL::Array)::Array
	id = Matrix{Float64}(I,length(hcen[1,:]),length(hcen[1,:]))
	inV = inv(hcen-1im*0.5*GammaT-chiLP*id)
	return etaLP/beta .* inV * GammaL
end


function Current(CurrentMatrix::Array)::Float64
	return real(tr(CurrentMatrix))
end


function DensityFExOpen(lx  ::Int,
                    ly  ::Int,
                    hcen::Array,
		    etaLP::Array,
		    chiLP::Array,
		    GammaT::Array,
		    beta ::Float64)::Array
	id = Matrix{Float64}(I,Int(lx*ly*2),Int(lx*ly*2))
	function DensMatInv(chilp,hcen,GammaT,etalp,beta)
		inV = inv(chilp*id-hcen+0.5*1im*GammaT)
		sum = real(etalp/beta * inV)
		return sum
	end
	rho = 1/2 * id + sum([DensMatInv(chiLP[1][i],hcen,GammaT,etaLP[i],beta) for i in 1:length(chiLP[1])]) + sum([DensMatInv(chiLP[2][i],hcen,GammaT,etaLP[i],beta) for i in 1:length(chiLP[2])])
    return rho
end

function SpinCurrent(CurrentMatrix::Array{ComplexF64})::Array{Float64}
	n = length(CurrentMatrix[:,1])
	a_sigmax = zeros(n,n)*(1+0*1im)
	a_sigmay = zeros(n,n)*(1+0*1im)
    	a_sigmaz = zeros(n,n)*(1+0*1im)
	for i in 1:2:n
		a_sigmax[i,i+1] = Sigma_x[1,2]
		a_sigmax[i+1,i] = Sigma_x[2,1]		
		a_sigmay[i,i+1] = Sigma_y[1,2]
		a_sigmay[i+1,i] = Sigma_y[2,1]
		a_sigmaz[i,i] 	= Sigma_z[1,1]
		a_sigmaz[i+1,i+1] = Sigma_z[2,2]
	end
	sx::Float64 = 0.0 
	sy::Float64 = 0.0
	sz::Float64 = 0.0
	cX = CurrentMatrix*a_sigmax
	cY = CurrentMatrix*a_sigmay
	cZ = CurrentMatrix*a_sigmaz
	for i in 1:n
		sx += real(cX[i,i])
		sy += real(cY[i,i])
		sz += real(cZ[i,i])
	end
		
	return [sx,sy,sz]
end


