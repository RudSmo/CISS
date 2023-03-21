# THIS IS USING THE AUXILIARY-MODE EXPANSION 
# BY POPESCU AND CROY (2016, New J. Phys. 18)

module Electrons
include("global_constants.jl")
include("pade.jl")
using .GlobalConstants
using .Pade
using LinearAlgebra

export OmegaEOM, PsiEOM, CurrentMatrix, RhoEOM, Rho0_Equilibrium_Closed, Rho0_Equilibrium_Open

function OmegaEOM(nAll,ChiPade,DeltaL,Omega0)
	nSys,nLeads,nPoles,nChannels = nAll
	Omega = zeros(nLeads,nPoles,nChannels,nLeads,nPoles,nChannels) * (1+0.0*1im)
	for a in 1:nLeads
		for p in 1:nPoles
			for c in 1:nChannels
				for ap in 1:nLeads
					for pp in 1:nPoles
						for cp in 1:nChannels
							Omega[a,p,c,ap,pp,cp] = -1im * (ChiPade[ap,pp,2] + DeltaL[ap] - ChiPade[a,p,1] - DeltaL[a]) * Omega0[a,p,c,ap,pp,cp]
						end
					end
				end
			end
		end
	end
	return Omega
end

function PsiEOM(nAll,H,DeltaL,BetaL,ChiPade,EtaPade,Gamma,Axi,Omega,Psi0)
	nSys,nLeads,nPoles,nChannels = nAll
	id = Matrix{Float64}(I,size(H))
	Apsi = zeros(nLeads,nPoles,nChannels,nSys) * (1+0.0*1im)
	
	for a in 1:nLeads
		for p in 1:nPoles
			A = H - (ChiPade[a,p,1] + DeltaL[a]) * id
			for c in 1:nChannels	
				B = EtaPade[p]/BetaL[a] * Gamma[a,c] * Axi[a,c,:] 
				C = zeros(nSys) * (1+0.0*1im)
				for ap in 1:nLeads
					for pp in 1:nPoles
						for cp in 1:nChannels
							C += Omega[a,p,c,ap,pp,cp] * Axi[ap,cp,:]
						end
					end
				end
				Apsi[a,p,c,:] = -1im * ( A * Psi0[a,p,c,:] + B + C )
			end
		end
	end
	return Apsi
end


function CurrentMatrix(nAll,Apsi,Axi)
	# Apsi = Array(nLeads,nPoles,nChannels,nSys)
	# Axi = Array(nLeads,nChannels,nSys)
	
	nSys,nLeads,nPoles,nChannels = nAll
	# CurrentMatrices for all leads			
	Api = zeros(nLeads,nSys,nSys)*(1+0.0*1im)     

	for a in 1:nLeads
		for c in 1:nChannels
			for p in 1:nPoles
				Api[a,:,:] += Apsi[a,p,c,:]*Axi[a,c,:]'
			end
		end
	end
	return Api
end


function RhoEOM(nAll,H,APi,rho)
	nSys,nLeads,nPoles,nChannels = nAll
	A = -1im * (H*rho-rho*H)
	B = zeros(size(rho)) * (1+0.0*1im)
	for a in 1:nLeads
		B[:,:] += APi[a,:,:] + conj(transpose(APi[a,:,:]))
	end
	return A + B
end




function Rho0_Equilibrium_Closed(L  ::Int,
                    		hcen::Array,
                    		fil :: Int,
                    		ex  :: Int)
        eigvec = eigvecs(hcen)
        teigvc = conj(transpose(eigvec))
        rho = zeros(2*L,2*L)*(1+1im*0)
        for i in 1:Int(2*L)
                for j in 1:Int(2*L)
                        rho[i,j] = sum([eigvec[i,a]*teigvc[a,j] for a in 1:(fil)])#+eigvec[i,fil+ex]*teigvc[fil+ex,j]
                end
        end
        return rho
end

function Rho0_Equilibrium_Open(L  	::Int,
                    	hcen	::Array,
			nPoles	::Int,
		    	EtaPade	::Array,
		    	ChiPade	::Array,
		    	GammaT	::Array,
		    	BetaL 	::Array)
	id = Matrix{Float64}(I,Int(2*L),Int(2*L))
	function DensMatInv(chilp,hcen,GammaT,etalp,beta)
		inV = inv(chilp*id-hcen+0.5*1im*GammaT)
		sum = real(etalp/beta * inV)
		return sum
	end
	rho = 1/2 * id + sum([DensMatInv(ChiPade[1,i,1],hcen,GammaT,EtaPade[i],BetaL[1]) for i in 1:nPoles]) + sum([DensMatInv(ChiPade[2,i,1],hcen,GammaT,EtaPade[i],BetaL[2]) for i in 1:nPoles])
    return rho
end




end
