module ClassicalPhonons
include("hamiltonian.jl")
using LinearAlgebra

using .Hamiltonian
export InitialPh,PhononPopulation,SampleWignerA,SymplecticEulerStep,SymplecticEulerStepOld,SymplecticEulerStepSingle

function InitialPh(popPh,ws,phase)
    return [sqrt(2*popPh[i]/ws[i])*cos(phase[i]) for i in 1:length(popPh)],
           [-sqrt(2*popPh[i]*ws[i])*sin(phase[i]) for i in 1:length(popPh)]
end
function PhononPopulation(qn,pn,ws)
    return [qn[i]^2*ws[i]/2 + pn[i]^2/(2*ws[i]) for i in 1:length(qn)]
end
function SampleWignerA(beta,w,version)
    """
    version == "v1": Distribution as given in eq. 51 of J. Chem. Phys. 152, 084110
    version == "v2": Distribution as given in eq. 14 of J. Chem. Phys. 156, 024108
    """
    aj = tanh(0.5*beta*w)
    if version == "v2"
        A  = randn(2)*aj/pi
    elseif version == "v1"
        A = randn(2)*2*aj
    end
    return [A[1]/sqrt(w*aj),A[2]*sqrt(w/aj)]
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

function SymplecticEulerStepOld(N,a,c,num_laps,tph,tsoph,ws,rho,q0,p0,dt)
    #A = rho * dHelph(N,a,c,num_laps,tph,tsoph,ws[1],1)
    #println(tr(A[1:2,1:2]))
    pn = real.([p0[i] - dt*(ws[i]^2 * q0[i] + tr(rho * dHelph(N,a,c,num_laps,tph,tsoph,ws[i],i))) for i in 1:N])
    qn = real.([q0[i] + dt*pn[i] for i in 1:N])
    return qn,pn
end
function SymplecticEulerStep(N,a,c,num_laps,tph,tsoph,ws,rho,q0,p0,dt)
    pn = []
    qn = []
    for i in 1:N
        A = rho * dHelph(N,a,c,num_laps,tph,tsoph,ws[i],i)
	pn_c = p0[i] - dt*(ws[i]^2 * q0[i] + tr(A[2*i-1:2*i,2*i-1:2*i]))
	qn_c = q0[i] + dt*pn_c
	push!(pn,real(pn_c))
	push!(qn,real(qn_c))
    end
    return qn,pn
end


function SymplecticEulerStepSingle(ws,gph,rho,q0,p0,dt)
    # println(tr(rho * dHelph(N,a,c,num_laps,tph,tsoph,ws[1],1)))
     pn = real(p0 - dt*(ws^2 * q0 + tr(rho * sqrt(2*ws)*gph*Matrix(I,2,2))))
     qn = real(q0 + dt*pn)
     return qn,pn
 end
end
