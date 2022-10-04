module NEGF
include("auxiliary.jl")
using .Auxiliary
using LinearAlgebra
using PyPlot
export CurSingle_A,LB_Transmission_A_single,GammaMatrix, G_R

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

function LB_Transmission_A_single(E,sys_par,lead_par,eta)            

    
    GammaLu,GammaLd,GammaLud,GammaLdu = SpinResolvedMatrix(GammaMatrix(sys_par,lead_par,"L"))
    GammaRu,GammaRd,GammaRud,GammaRdu = SpinResolvedMatrix(GammaMatrix(sys_par,lead_par,"R")) 
    SigmaR                            = -1im/2 .* (GammaMatrix(sys_par,lead_par,"L")+GammaMatrix(sys_par,lead_par,"R"))
    GR                                = G_R(sys_par.HS,E,eta,SigmaR)
    GRu,GRd,GRud,GRdu                 = SpinResolvedMatrix(GR)
    
    return LB_Transmission(GRu,GammaLu,GammaRu), 
           LB_Transmission(GRd,GammaLd,GammaRd), 
           LB_Transmission(GRud,GammaLu,GammaRd), 
           LB_Transmission(GRdu,GammaLd,GammaRu)
end


function CurSingle_A(i,sys_par,lead_par,eta)
    E1,E2 = sys_par.En[i:i+1]
    dE    = E2-E1
    T1u,T1d,T1ud,T1du = LB_Transmission_A_single(E1,sys_par,lead_par,eta) 
    T2u,T2d,T2ud,T2du = LB_Transmission_A_single(E2,sys_par,lead_par,eta) 
   
    fermi1u  = Fermi(E1,lead_par.muLu,lead_par.tLu) - Fermi(E1, lead_par.muRu,lead_par.tRu)
    fermi1d  = Fermi(E1,lead_par.muLd,lead_par.tLd) - Fermi(E1, lead_par.muRd,lead_par.tRd)
    fermi1ud = Fermi(E1,lead_par.muLd,lead_par.tLd) - Fermi(E1, lead_par.muRu,lead_par.tRu)
    fermi1du = Fermi(E1,lead_par.muLu,lead_par.tLu) - Fermi(E1, lead_par.muRd,lead_par.tRd)
    
    fermi2u  = Fermi(E2,lead_par.muLu,lead_par.tLu) - Fermi(E2, lead_par.muRu,lead_par.tRu)
    fermi2d  = Fermi(E2,lead_par.muLd,lead_par.tLd) - Fermi(E2, lead_par.muRd,lead_par.tRd)
    fermi2ud = Fermi(E2,lead_par.muLd,lead_par.tLd) - Fermi(E2, lead_par.muRu,lead_par.tRu)
    fermi2du = Fermi(E2,lead_par.muLu,lead_par.tLu) - Fermi(E2, lead_par.muRd,lead_par.tRd)


    f1u = T1u * fermi1u
    #println("$(T1u) $(fermi1u) $(f1u)\n")
    f2u = T2u * fermi2u
    
    f1d = T1d * fermi1d
    f2d = T2d * fermi2d 

    f1ud = T1ud*fermi1ud
    f2ud = T2ud*fermi2ud 

    f1du = T1du*fermi1du
    f2du = T2du*fermi2du 
    return dE*0.5*(f1u+f2u),dE*0.5*(f1d+f2d),dE*0.5*(f1ud+f2ud),dE*0.5*(f1du+f2du)
end






end