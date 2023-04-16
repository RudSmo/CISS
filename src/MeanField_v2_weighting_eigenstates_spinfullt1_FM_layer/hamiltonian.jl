module Hamiltonian
using LinearAlgebra

export VDot,absv,UnitVector,FranssonHamiltonian,FranssonHamiltonian_H0,FranssonHamiltonian_HI,FranssonHamiltonianNH,HolsteinHopping,SpinHopping


function VDot(a)
    sx = [0 1; 1 0.]
    sy = [0 -1im; 1im 0.]
    sz = [1 0; 0 -1.] 
    return a[1]*sx + a[2]*sy + a[3]*sz
end
function absv(a)
    return sqrt(dot(a,a))
end
#-------------------------- Electronic Hamiltonian-------------------------------------
function UnitVector(a,c,j,num_laps,N)
    phij = 2*pi*(j-1)*num_laps/(N-1)
    cj = c*(j-1)/(N-1)
    return [a*cos(phij),a*sin(phij),cj]
end
function SOC_A_nn(N,mol_sites,tso,a,c,num_laps)
    Hs = zeros(2*N,2*N)*(1+0.0*1im)
    for k in mol_sites[1:end-2]
	j = 2*k-1
        rj = UnitVector(a,c,k,num_laps,N)
        rj1 = UnitVector(a,c,k+1,num_laps,N)
        rj2 = UnitVector(a,c,k+2,num_laps,N) 
        rjrj1 = rj-rj1
        rjrj2 = rj-rj2

        dj1 = rjrj1/absv(rjrj1)
        dj2 = rjrj2/absv(rjrj2)
        
        vj = cross(dj1,dj2)
         
        Vj = 1im*tso*VDot(vj)

        Hs[j,j+4] = Vj[1,1]
        Hs[j+1,j+4] = Vj[2,1]
        Hs[j,j+5] = Vj[1,2]
        Hs[j+1,j+5] = Vj[2,2]  
    end
    return Hs+Hs'
end

function Onsite_A(N,subSites,E)
    eu,ed = E
    H = zeros(2*N,2*N)
    for i in subSites
	j = 2*i-1
        H[j,j] = eu
        H[j+1,j+1] = ed
    end
    return H
end

function Hopping_A_nn(N,subSites,t)
    Hh = zeros(2*N,2*N)
    for i in subSites[1:end-1]
	Hh[i,i+2] = t
       	Hh[i+2,i] = t
    end
    return Hh
end

function FM_Kondo(N,fm_sites,Jsd,mag)
    H = zeros(2*N,2*N)
    j = 1
    for i in fm_sites
	hm = Jsd * Vdot(mag[j])
	H[i,i] = hm[1,1]
	H[i,i+1] = hm[1,2]
	H[i+1,i] = hm[2,1]
	H[i+1,i+1] = hm[2,2]
	j += 1
    end
    return H
end

function FM_molecule_coupl(N,coup_site,coup_en)
	H = zeros(2*N,2*N)
	i = 2*coup_site-1
	H[i,i+2] = coup_en[1,1]
	H[i,i+3] = coup_en[1,2]	
	H[i+1,i+2] = coup_en[2,1]
	H[i+1,i+3] = coup_en[2,2]
	
	H = H + H'
	return H
end

#-------------------------- Electron-Phonon Hamiltonian-------------------------------------

function Peierls(N,mol_sites,tph,ws,xs)
    H = zeros(2*N,2*N)
    for k in mol_sites[1:end-1]
	i = 2*k-1
        H[i,i+2] = tph[1] * sqrt(2*ws[k])*xs[k] 
	H[i+1,i+3] = tph[2] * sqrt(2*ws[k])*xs[k] 
        H[i+2,i] = tph[1] * sqrt(2*ws[k])*xs[k]  
	H[i+3,i+1] = tph[2] * sqrt(2*ws[k])*xs[k]
    end
    return H
end


function Peierls_SOC(N,mol_sites,tsoph,a,c,num_laps,ws,xs)
    Hs = zeros(2*N,2*N)*(1+0.0*1im)
    for k in mol_sites[1:end-2]
	j = 2*k-1
        rj = UnitVector(a,c,k,num_laps,N)
        rj1 = UnitVector(a,c,k+1,num_laps,N)
        rj2 = UnitVector(a,c,k+2,num_laps,N) 
        rjrj1 = rj-rj1
        rjrj2 = rj-rj2

        dj1 = rjrj1/absv(rjrj1)
        dj2 = rjrj2/absv(rjrj2)
        
        vj = cross(dj1,dj2) 
        Vj = 1im*tsoph*VDot(vj) * sqrt(2*ws[k]) * xs[k]

        Hs[j,j+4] = Vj[1,1]
        Hs[j+1,j+4] = Vj[2,1]
        Hs[j,j+5] = Vj[1,2]
        Hs[j+1,j+5] = Vj[2,2]  
    end
    return Hs+Hs'
end



function FranssonHamiltonianFM(HamPar)
	# FM layer
        HoFM = Onsite_A(HamPar.N,HamPar.fm_sites,HamPar.OnsiteEnergiesFM) 
        HhFM = Hopping_A_nn(HamPar.N,HamPar.fm_sites,HamPar.Hopping_FM)
	HKFM = FM_Kondo(HamPar.N,HamPar.fm_sites,HamPar.Jsd,HamPar.mag)
	
	# FM-molecule coupling
	HcFM_mol = FM_molecule_coupl(HamPar.N,HamPar.coup_site,HamPar.coup_en)
	
	# Molecule	
	Ho = Onsite_A(HamPar.N,HamPar.mol_sites,HamPar.OnsiteEnergies)
	Hh = Hopping_A_nn(HamPar.N,HamPar.mol_sites,HamPar.Hopping)
        HSo = SOC_A_nn(HamPar.N,HamPar.mol_sites,HamPar.SOC,HamPar.a_hel,HamPar.c_hel,HamPar.num_laps)
        Heph = Peierls(HamPar.N,HamPar.mol_sites,HamPar.HoppingPh,HamPar.ws,HamPar.X)
        Hsoph = Peierls_SOC(HamPar.N,HamPar.mol_sites,HamPar.SOCPh,HamPar.a_hel,HamPar.c_hel,HamPar.num_laps,HamPar.ws,HamPar.X)
        return HoFM + HhFM + HKFM + HcFM_mol + Ho+Hh + HSo + Heph + Hsoph
end

end
