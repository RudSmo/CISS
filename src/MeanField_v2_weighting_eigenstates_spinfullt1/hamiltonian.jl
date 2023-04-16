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
function SOC_A_nn(N,tso,a,c,num_laps)
    Hs = zeros(2*N,2*N)*(1+0.0*1im)
    for k in 1:N-2
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

function Hopping_A_nn(N,t)
    Hh = zeros(2*N,2*N)
    for i in 1:2*N-2
        Hh[i,i+2] = t
        Hh[i+2,i] = t
    end
    return Hh
end
function Hopping_A_nn_spin(N,t)
    tu,td = t
    Hh = zeros(2*N,2*N)
    for i in 1:2*N-2
        Hh[i,i+2] = tu
        Hh[i+2,i] = td
    end
    return Hh
end
function Onsite_A(N,E)
    eu,ed = E 
    H = zeros(2*N,2*N)
    for i in 1:2:2*N-1
        H[i,i] = eu 
    end
    for i in 2:2:2*N
        H[i,i] = ed 
    end
    return H
end


#-------------------------- Electron-Phonon Hamiltonian-------------------------------------
function Holstein(N,gPh,ws,xs)
    H = zeros(2*N,2*N)
    for k in 1:N
    	i = 2*k-1    
	H[i,i] = gPh[1] * sqrt(2*ws[k])*xs[k]
        H[i+1,i+1] = gPh[2] * sqrt(2*ws[k])*xs[k] 
    end
    return H
end


function Peierls_Old(N,tph,ws,xs)
    H = zeros(2*N,2*N)
    for k in 1:N-1
	i = 2*k-1
        H[i,i+2] = tph[1] * sqrt(2*ws[k])*xs[k] 
	H[i+1,i+3] = tph[2] * sqrt(2*ws[k])*xs[k] 
        H[i+2,i] = tph[1] * sqrt(2*ws[k])*xs[k]  
	H[i+3,i+1] = tph[2] * sqrt(2*ws[k])*xs[k]
    end
    return H
end


function Peierls(N,tph,ws,xs)
    H = zeros(2*N,2*N)
    k = 1
    for i in 1:2*N-2
        H[i,i+2] = tph * sqrt(2*ws[k])*xs[k]  
        H[i+2,i] = tph * sqrt(2*ws[k])*xs[k]  
        if i % 2 == 0
            k += 1
        end
    end
    H[2*N,2*N-2] = tph * sqrt(2*ws[N])*xs[N]
    H[2*N-2,2*N] = tph * sqrt(2*ws[N])*xs[N]
    return H
end



function Peierls_SOC(N,tsoph,a,c,num_laps,ws,xs)
    Hs = zeros(2*N,2*N)*(1+0.0*1im)
    for k in 1:N-2
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



function FranssonHamiltonian(HamPar)
	Ho = Onsite_A(HamPar.N,HamPar.OnsiteEnergies)
	Hh = Hopping_A_nn(HamPar.N,HamPar.Hopping)
	HSo = SOC_A_nn(HamPar.N,HamPar.SOC,HamPar.a_hel,HamPar.c_hel,HamPar.num_laps)
	#if ver == "0"
    	Heph = Peierls_Old(HamPar.N,HamPar.HoppingPh,HamPar.ws,HamPar.X)
    	#Heph = Peierls(HamPar.N,HamPar.HoppingPh,HamPar.ws,HamPar.X)
	#else
	#	Heph = Peierls(N,tph,ws,xs)
	#end
    	Hsoph = Peierls_SOC(HamPar.N,HamPar.SOCPh,HamPar.a_hel,HamPar.c_hel,HamPar.num_laps,HamPar.ws,HamPar.X)
	return Ho+Hh+HSo+Heph+Hsoph
end
function FranssonHamiltonian_H0(HamPar)
        Ho = Onsite_A(HamPar.N,HamPar.OnsiteEnergies)
        Hh = Hopping_A_nn(HamPar.N,HamPar.Hopping)
        HSo = SOC_A_nn(HamPar.N,HamPar.SOC,HamPar.a_hel,HamPar.c_hel,HamPar.num_laps)
        return Ho+Hh+HSo
end

function FranssonHamiltonian_HI(HamPar)
    	Heph = Peierls_Old(HamPar.N,HamPar.HoppingPh,HamPar.ws,HamPar.X)
    	Hsoph = Peierls_SOC(HamPar.N,HamPar.SOCPh,HamPar.a_hel,HamPar.c_hel,HamPar.num_laps,HamPar.ws,HamPar.X)
        return Heph+Hsoph
end

function FranssonHamiltonianNH(HamPar)
    HF = FranssonHamiltonian(HamPar)
    HF[1,1] -= 1im/2* HamPar.GL[1]
    HF[2,2] -= 1im/2* HamPar.GL[2]
    HF[end-1,end-1] -= 1im/2* HamPar.GR[1]
    HF[end,end] -= 1im/2 * HamPar.GR[2]
    return HF
end

function HolsteinHopping(HamPar)
    Ho = Onsite_A(HamPar.N,HamPar.OnsiteEnergies)
    Hh = Hopping_A_nn(HamPar.N,HamPar.Hopping)
    Hph = Holstein(HamPar.N,HamPar.gPh,HamPar.ws,HamPar.X)
    Htph = Peierls_Old(HamPar.N,HamPar.HoppingPh,HamPar.ws,HamPar.X)
    return Ho .+ Hh .+ Hph .+ Htph
end

function SpinHopping(HamPar,t0)
	Ho = Onsite_A(HamPar.N,HamPar.OnsiteEnergies)
	Hh = Hopping_A_nn_spin(HamPar.N,t0)
	return Ho .+ Hh	
end

end
