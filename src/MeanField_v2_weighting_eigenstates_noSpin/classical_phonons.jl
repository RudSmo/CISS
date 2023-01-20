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
	Pn = []
	Xn =  []

	for i in 1:HamPar.N
	    nd = rho[i,i] 
	    if i < HamPar.N
		    t_elmt = HamPar.HoppingPh * (rho[i,i+1] + rho[i+1,i])
	    else
		    t_elmt = 0.0	
	    end

	    Pn_c = HamPar.P[i] - dt * (sqrt(2*HamPar.ws[i]) * (HamPar.gPh *  nd + t_elmt) + HamPar.ws[i]^2 * HamPar.X[i])
	    Xn_c = HamPar.X[i] + dt * Pn_c
	    push!(Pn,real(Pn_c))
	    push!(Xn,real(Xn_c))
	end
        return Xn,Pn
end

end
