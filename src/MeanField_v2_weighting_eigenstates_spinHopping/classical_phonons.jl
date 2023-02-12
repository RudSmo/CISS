module ClassicalPhonons
include("hamiltonian.jl")
using LinearAlgebra
using PyCall
random = pyimport("random")
using .Hamiltonian
export InitialPh,PhononPopulation,SampleWigner,SampleWigner_Multi,Wignerweights,SymplecticEulerStep

function InitialPh(popPh,ws,phase)
    return [sqrt(2*popPh[i]/ws[i])*cos(phase[i]) for i in 1:length(popPh)],
           [-sqrt(2*popPh[i]*ws[i])*sin(phase[i]) for i in 1:length(popPh)]
end

function PhononPopulation(qn,pn,ws)
    return [qn[i]^2*ws[i]/2 + pn[i]^2 for i in 1:length(qn)]
end

function Wigner(X,P,beta,w)
    W = zeros(length(P),length(X))
    for i in 1:length(P)
        for j in 1:length(X)
            W[i,j] = 1/pi * tanh(0.5*beta*w) * exp(-tanh(0.5*beta*w)/w * (w^2 * X[j]^2 + P[i]^2) )
        end
    end
    return W
 end

function SampleWigner(beta,w,version)
    """
    version == "v1": Distribution as given in eq. 51 of J. Chem. Phys. 152, 084110
    version == "v2": Distribution as given in eq. 14 of J. Chem. Phys. 156, 024108
    """
    aj = tanh(0.5*beta*w)

    if version == "v1"
        A = randn(2)*2*aj
        return [A[1]/sqrt(w*aj),A[2]*sqrt(w/aj)]
    elseif version == "v2"
        A  = randn(2)*aj/pi
        return [A[1]/sqrt(w*aj),A[2]*sqrt(w/aj)]
    elseif version == "new"
        A = randn(2)
        return [A[1]/sqrt(w*aj),A[2]*sqrt(w/aj)] 
    elseif version == "correct"
        # This implementation is really slow, because the analytical Wigner probability density has to be generated for each vibrational mode
        P = Array(-1:0.002:1)
        X = Array(-15:0.01:15)
        W = Wigner(X,P,beta,w)
        c = [(random.choices(P,W[:,j],k=m),random.choices(X,W[i,:],k=m)) for i in [random.randint(1,length(P)) for k in 1:m] for j in [random.randint(1,length(X)) for k in 1:m]];
        return [c[1][2],c[1][1]]
    end
end

function SampleWigner_Multi(beta,w,m)
    # Assume w[i] = w[1]
    P = Array(-1:0.002:1)
    X = Array(-15:0.02:15) 
    aj = tanh(0.5*beta*w[1])
    W = Wigner(X,P,beta,w[1]) 
    Xr = []
    Pr = []
    for i in 1:m
        push!(Xr,random.choices(X,W[random.randint(1,length(P)),:],k=1)[1])
        push!(Pr,random.choices(P,W[:,random.randint(1,length(X))],k=1)[1])
    end
    return Xr,Pr
end



function Wignerweights(X,P,beta,ws,version)
    if version == "v2" || version == "new"
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
	    if i < HamPar.N
	            t_elmt = HamPar.HoppingPh* tr(rho[2*i-1:2*i,2*i+1:2*i+2] + rho[2*i+1:2*i+2,2*i-1:2*i])
	    else
		    t_elmt = 0.0
	    end
            Pn_c = HamPar.P[i] - dt * (sqrt(2*HamPar.ws[i]) * (HamPar.gPh[1] *  nd_up + HamPar.gPh[2] *  nd_dn + t_elmt) + HamPar.ws[i]^2 * HamPar.X[i])
	    Xn_c = HamPar.X[i] + dt * Pn_c
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
	        Xn_c = HamPar.X[i] + dt * Pn_c
            elseif i <= HamPar.N-1
                t_elmt = HamPar.HoppingPh*tr(rho[2*i-1:2*i,2*i+1:2*i+2] + rho[2*i+1:2*i+2,2*i-1:2*i])
                Pn_c = HamPar.P[i] - dt * (sqrt(2*HamPar.ws[i]) * t_elmt + HamPar.ws[i]^2 * HamPar.X[i])
	        Xn_c = HamPar.X[i] + dt * Pn_c 
            else
                Pn_c = HamPar.P[i] - dt * HamPar.ws[i]^2 * HamPar.X[i]
                Xn_c = HamPar.X[i] + dt * Pn_c
            end 
	        push!(Pn,real(Pn_c))
	        push!(Xn,real(Xn_c))
        end 
        return Xn,Pn
    end
end

end
