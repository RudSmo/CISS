module NEGF
include("auxiliary.jl")
using .Auxiliary
using LinearAlgebra
export CurSingle,LB_Transmission,GammaMatrix, G_R

function GammaMatrix(sys_par,lead_par,lead)
    Ga = zeros(size(sys_par.HS))*(1+0.0*1im)
    if lead == "L"
        Ga[1,1] = lead_par.gLu
        Ga[2,2] = lead_par.gLd
        Ga[1,2] = lead_par.gLud
        Ga[2,1] = lead_par.gLdu
    else
        Ga[end-1,end-1] = lead_par.gRu
        Ga[end,end]     = lead_par.gRd 
        Ga[end-1,end]   = lead_par.gRud
        Ga[end,end-1]   = lead_par.gRdu
    end
    return Ga
end

function G_R(HS,E,eta,Sigma_R)
    id = Matrix(I,size(HS))
    return inv(E * id - 1im * eta * id - HS - Sigma_R)
end

function LB_Transmission(GR,GammaL,GammaR)
    GA = GR'
    return 1/(2*pi) * tr(GammaL*GR*GammaR*GA)
end

function CurSingle(i,sys_par,lead_par,eta)
    E1,E2 = sys_par.En[i:i+1]
    GammaL= GammaMatrix(sys_par,lead_par,"L")
    GammaR= GammaMatrix(sys_par,lead_par,"R")
    Sigma_R = -1im/2 .* (GammaL+GammaR) 
    dE    = E2-E1
    HS    = sys_par.HS
    GR1   = G_R(HS,E1,eta,Sigma_R)
    GR2   = G_R(HS,E2,eta,Sigma_R)
    T1    = LB_Transmission(GR1,GammaL,GammaR) 
    T2    = LB_Transmission(GR2,GammaL,GammaR) 
   
    fermi1  = Fermi(E1,lead_par.muL,lead_par.tL) - Fermi(E1, lead_par.muR,lead_par.tR)    
    fermi2  = Fermi(E2,lead_par.muL,lead_par.tL) - Fermi(E2, lead_par.muR,lead_par.tR)


    f1 = T1 * fermi1
    f2 = T2 * fermi2
    
    return dE*0.5*(f1+f2)
end






end
