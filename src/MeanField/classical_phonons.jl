module ClassicalPhonons
include("hamiltonian.jl")
using LinearAlgebra

using .Hamiltonian
export SampleWignerA,SymplecticEulerStep

function SampleWignerA(beta,w)
    A= randn(2)*w*beta
    return [A[1]/w,A[2]]
end

function dHelph(N,a,c,num_laps,tph,tsoph,w,n)
    dH = zeros(2*N,2*N)*(1+0.0*1im)
    k = Int(2*n-1)

    if n <= N-2
        rj = UnitVector(a,c,n,num_laps,N)
        rj1 = UnitVector(a,c,n+1,num_laps,N)
        rj2 = UnitVector(a,c,n+2,num_laps,N) 
        rjrj1 = rj-rj1
        rjrj2 = rj-rj2

        dj1 = rjrj1/absv(rjrj1)
        dj2 = rjrj2/absv(rjrj2)
    
        vj = cross(dj1,dj2) 
        Vj = 1im*tsoph*VDot(vj) * sqrt(2*w) 

        dH[k,k+4] = Vj[1,1]
        dH[k+1,k+4] = Vj[2,1]
        dH[k,k+5] = Vj[1,2]
        dH[k+1,k+5] = Vj[2,2]  
    end
    if n <= N-1
        dH[k,k+2] = -tph*sqrt(2*w)
        dH[k+2,k] = -tph*sqrt(2*w)
        dH[k+1,k+3] = -tph*sqrt(2*w)
        dH[k+3,k+1] = -tph*sqrt(2*w)
    end
     
    return dH
end

function SymplecticEulerStep(N,a,c,num_laps,tph,tsoph,ws,rho,q0,p0,dt)
   # println(tr(rho * dHelph(N,a,c,num_laps,tph,tsoph,ws[1],1)))
    pn = real.([p0[i] - dt*(ws[i]^2 * q0[i] - tr(rho * dHelph(N,a,c,num_laps,tph,tsoph,ws[i],i))) for i in 1:N])
    qn = real.([q0[i] + dt*pn[i] for i in 1:N])
    return qn,pn
end
end