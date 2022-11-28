module ClassicalPhonons
include("hamiltonian.jl")
using LinearAlgebra

using .Hamiltonian
export InitialPh,PhononPopulation,SampleWigner,Wignerweights,SymplecticEulerStep

function InitialPh(popPh,ws,phase)
    return [sqrt(2*popPh[i]/ws[i])*cos(phase[i]) for i in 1:length(popPh)],
           [-sqrt(2*popPh[i]*ws[i])*sin(phase[i]) for i in 1:length(popPh)]
end

function PhononPopulation(qn,pn,ws)
    return [qn[i]^2*ws[i]/2 + pn[i]^2 for i in 1:length(qn)]
end

function SampleWigner(beta,w,version)
    """
    version == "v1": Distribution as given in eq. 51 of J. Chem. Phys. 152, 084110
    version == "v2": Distribution as given in eq. 14 of J. Chem. Phys. 156, 024108
    """
    aj = tanh(0.5*beta*w)
    if version == "v2"
        A  = randn(2)*aj/pi
    elseif version == "v1"
        A = randn(2)*2*aj
    elseif version == "new"
        A = randn(2)
    end
    return [A[1]/sqrt(w*aj),A[2]*sqrt(w/aj)]
end

function Wignerweights(X,P,beta,ws,version)
    if version == "v2"
        return [1/pi * tanh(0.5*beta*ws[i]) * exp(-tanh(0.5*beta*ws[i]) * (ws[i]* X[i]^2 + P[i]^2/ws[i])) for i in 1:length(X)]
    elseif version == "v1"
        return [2*tanh(0.5*beta*ws[i]) * exp(-tanh(0.5*beta*ws[i]) * (ws[i]* X[i]^2 + P[i]^2/ws[i])) for i in 1:length(X)]
    end
end

function SymplecticEulerStep(HamPar,rho,dt)
    if HamPar.HamiltonStr == "Holstein"
        Pn = []
        Xn =  []
    
        for i in 1:HamPar.N
            nd_up = rho[2*i-1,2*i-1] 
            nd_dn = rho[2*i,2*i] 
            Pn_c = HamPar.P[i] - dt * (sqrt(2*HamPar.ws[i]) * (HamPar.gPh[1] *  nd_up + HamPar.gPh[2] *  nd_dn) + HamPar.ws[i]^2 * HamPar.X[i])
	        Xn_c = HamPar.X[i] + dt * HamPar.ws[i] * Pn_c
	        push!(Pn,real(Pn_c))
	        push!(Xn,real(Xn_c))
        end
        return Xn,Pn
    elseif HamPar.HamiltonStr == "Fransson"
        Pn = []
        Xn = []
        for i in 1:HamPar.N
            if i <= HamPar.N-2
                rj = UnitVector(HamPar.a_hel,HamPar.c_hel,i,HamPar.num_laps,HamPar.N)
                rj1 = UnitVector(HamPar.a_hel,HamPar.c_hel,i+1,HamPar.num_laps,HamPar.N)
                rj2 = UnitVector(HamPar.a_hel,HamPar.c_hel,i+2,HamPar.num_laps,HamPar.N) 
                rjrj1 = rj-rj1
                rjrj2 = rj-rj2

                dj1 = rjrj1/absv(rjrj1)
                dj2 = rjrj2/absv(rjrj2)
        
                vj = cross(dj1,dj2) 
                coupl = 1im * HamPar.SOCPh

                so_elmt = coupl*tr((rho[2*i-1:2*i,2*i+3:2*i+4] + rho[2*i+3:2*i+4,2*i-1:2*i])  * [1 1-1im;1+1im -1] )
                t_elmt = HamPar.HoppingPh* tr(rho[2*i-1:2*i,2*i+1:2*i+2] + rho[2*i+1:2*i+2,2*i-1:2*i]) 
                Pn_c = HamPar.P[i] - dt * (sqrt(2*HamPar.ws[i]) * (t_elmt + so_elmt) + HamPar.ws[i]^2 * HamPar.X[i])
	            Xn_c = HamPar.X[i] + dt * HamPar.ws[i] * Pn_c
            elseif i <= HamPar.N-1
                t_elmt = HamPar.HoppingPh*tr(rho[2*i-1:2*i,2*i+1:2*i+2] + rho[2*i+1:2*i+2,2*i-1:2*i])
                Pn_c = HamPar.P[i] - dt * (sqrt(2*HamPar.ws[i]) * t_elmt + HamPar.ws[i]^2 * HamPar.X[i])
	            Xn_c = HamPar.X[i] + dt * HamPar.ws[i] * Pn_c 
            else
                Pn_c = HamPar.P[i] - dt * HamPar.ws[i]^2 * HamPar.X[i]
                Xn_c = HamPar.X[i] + dt * HamPar.ws[i] * Pn_c
            end 
	        push!(Pn,real(Pn_c))
	        push!(Xn,real(Xn_c))
        end 
        return Xn,Pn
    end
end

end
