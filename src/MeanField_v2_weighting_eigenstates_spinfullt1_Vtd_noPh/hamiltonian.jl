module Hamiltonian
using LinearAlgebra

export VDot,absv,UnitVector,FranssonHamiltonian


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
    k = 1
    for j in 1:2:2*N-4
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
        k += 1
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


function FranssonHamiltonian(HamPar)
        Ho = Onsite_A(HamPar.N,HamPar.OnsiteEnergies)
        Hh = Hopping_A_nn(HamPar.N,HamPar.Hopping)
        HSo = SOC_A_nn(HamPar.N,HamPar.SOC,HamPar.a_hel,HamPar.c_hel,HamPar.num_laps)
        return Ho+Hh+HSo
end
end
