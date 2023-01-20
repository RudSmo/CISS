module Hamiltonian
using LinearAlgebra

export HolsteinHopping


function Hopping_A_nn(N,t)
    Hh = zeros(N,N)
    for i in 1:N-1
        Hh[i,i+1] = -t
        Hh[i+1,i] = -t
    end
    return Hh
end
function Onsite_A(N,E) 
    H = zeros(N,N)
    for i in 1:N
        H[i,i] = E 
    end
    return H
end


#-------------------------- Electron-Phonon Hamiltonian-------------------------------------
function Holstein(N,gPh,ws,xs)
    H = zeros(N,N)
    k = 1
    for i in 1:N
        H[i,i] = gPh * sqrt(2*ws[k])*xs[k]
    end
    return H
end


function Peierls_Old(N,tph,ws,xs)
    H = zeros(N,N)
    k = 1
    for i in 1:N-1
        H[i,i+1] = -tph * sqrt(2*ws[k])*xs[k]  
        H[i+1,i] = -tph * sqrt(2*ws[k])*xs[k]  
    end
    return H
end


function Peierls(N,tph,ws,xs)
    H = zeros(N,N)
    k = 1
    for i in 1:N-1
        H[i,i+1] = -tph * sqrt(2*ws[k])*xs[k]  
        H[i+1,i] = -tph * sqrt(2*ws[k])*xs[k]  
    end
    H[N,N-1] = -tph * sqrt(2*ws[N])*xs[N]
    H[N-1,N] = -tph * sqrt(2*ws[N])*xs[N]
    return H
end


function HolsteinHopping(HamPar)
    Ho = Onsite_A(HamPar.N,HamPar.OnsiteEnergies)
    Hh = Hopping_A_nn(HamPar.N,HamPar.Hopping)
    Hph = Holstein(HamPar.N,HamPar.gPh,HamPar.ws,HamPar.X)
    Htph = Peierls_Old(HamPar.N,HamPar.HoppingPh,HamPar.ws,HamPar.X)
    return Ho .+ Hh .+ Hph .+ Htph
end

end
