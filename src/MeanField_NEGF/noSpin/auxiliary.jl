module Auxiliary
using LinearAlgebra
export SpinResolvedMatrix,Fermi

function SpinResolvedMatrix(H)
    sH = Int(size(H)[1]/2)
    Hu = zeros(sH,sH)*1im
    Hd = zeros(sH,sH)*1im
    Hud = zeros(sH,sH)*1im
    Hdu = zeros(sH,sH)*1im

    k = 1
    for i in 1:sH
        l = 1
        for j in 1:sH
            Hu[i,j] = H[k,l]
            l += 2
        end
        k += 2
    end

    k = 2
    for i in 1:sH
        l = 2
        for j in 1:sH
            Hd[i,j] = H[k,l]
            l += 2
        end
        k += 2
    end

    k = 1
    for i in 1:sH
        l = 2
        for j in 1:sH
            Hud[i,j] = H[k,l]
            l += 2
        end
        k += 2
    end

    k = 2
    for i in 1:sH
        l = 1
        for j in 1:sH
            Hdu[i,j] = H[k,l]
            l += 2
        end
        k += 2
    end

    return Hu,Hd,Hud,Hdu            
end

function Fermi(E,mu,T)
    return 1/(exp((E-mu)/T)+1)
end

function SemiElliptic(En,B,mu)
    return [1/(pi*B^2)*sqrt(B^2+(E-mu)) for E in En]
end

function SingleLorentzian(E,gamma,mu)
    return 1/((E^2-mu^2)^2+gamma^2*mu^2)
end

function Lorentzian(En,gamma,mu)
    if length(gamma) == length(mu)
        return [sum([SingleLorentzian(E,gamma[i],mu[i]) for i in 1:length(gamma)]) for E in En]
    else
        ArgumentError("length of gamma != length of mu. Provide as many energy-offsets as broadening parameters to the Lorentzian.")
    end
end



end