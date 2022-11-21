module Observables
include("negf.jl")
include("auxiliary.jl")
using .NEGF
using .Auxiliary
using LinearAlgebra
export SampleWigner,DOS, LB_CurrentV, LB_TransmissionE, LB_LocalCurrent, LB_LocalCurrentV

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
    end
    return [A[1]/sqrt(w*aj),A[2]*sqrt(w/aj)]
end
function DOS(En,HS,sigma)
    en = eigvals(HS)
    dos = []
    for i in 1:length(En) 
        ds = sum([1/(sigma*sqrt(pi))*exp(-(En[i]-ei)^2/(sigma^2)) for ei in en])
        push!(dos,ds)
    end
    return dos
end

function LB_TransmissionE(sys_par,lead_par,eta=1e-12)
    T = []
    GammaL = GammaMatrix(sys_par,lead_par,"L")
    GammaR = GammaMatrix(sys_par,lead_par,"R") 
    SigmaR = -1im/2 .* (GammaL+GammaR) 
    for E in sys_par.En
	GR = G_R(sys_par.HS,E,eta,SigmaR)
        t = LB_Transmission(GR,GammaL,GammaR)
        push!(T,t)
    end
    return T
end
function SPT(sys_par,lead_par)
    Tu,Td,Tud,Tdu = LB_TransmissionE(sys_par,lead_par) 
    return [(Tu[i]+Tdu[i]-Tud[i]-Td[i])/(Tu[i]+Tdu[i]+Tud[i]+Td[i]) for i in 1:length(Tu)]
end

function LB_Current(sys_par,lead_par,eta)
    """
    In model A the odd indices denote the up spin subspace, whereas the even ones the down spin subspace.
    """
    I = 0.0
    for i in 1:length(sys_par.En)-1
        ic = CurSingle(i,sys_par,lead_par,eta)
        I += real(ic)
    end
    return I
end

function LB_CurrentV(sys_par,lead_par,V,path,savestp,eta=1e-9)
    I = []
    mode= "w"
    k   = 0
    for (ind,v) in enumerate(V)
        lead_par.muL =  v/2
        lead_par.muR = - v/2
        ic = LB_Current(sys_par,lead_par,eta)
        push!(I,ic)
        if ind % savestp == 0
            fil = open(path,mode)
            if k == 0
                write(fil, "# V IL \n")
            end
            for (v,c) in zip(V[1+k*savestp:(k+1)*savestp],I)
                write(fil,"$v $ic \n")
            end
            I  = []
            k   += 1
	    mode = "a"
            close(fil)
        end
    end
end

function LB_LocalCurrent(sys_par,lead_par,eta=1e-9)
    Inu = zeros(sys_par.N-1) 
    Ind = zeros(sys_par.N-1)
    Inud= zeros(sys_par.N-1)
    Indu= zeros(sys_par.N-1)
    Hu,Hd,Hud,Hdu = SpinResolvedMatrix(sys_par.HS)     
    GammaL        = GammaMatrix(sys_par,lead_par,"L")
    GammaR        = GammaMatrix(sys_par,lead_par,"R")   
    GLu,GLd,GLud,GLdu = SpinResolvedMatrix(GammaL)
    GRu,GRd,GRud,GRdu = SpinResolvedMatrix(GammaR)
    SR                = -1im/2 * (GammaL + GammaR)
    for i in 1:length(sys_par.En)-1
        E1 = sys_par.En[i]
        E2 = sys_par.En[i+1] 
        dE = E2-E1
        fermi1u  = Fermi(E1,lead_par.muLu,lead_par.tLu) - Fermi(E1, lead_par.muRu,lead_par.tRu)    
        fermi1d  = Fermi(E1,lead_par.muLd,lead_par.tLd) - Fermi(E1, lead_par.muRd,lead_par.tRd)
        
        fermi2u  = Fermi(E2,lead_par.muLu,lead_par.tLu) - Fermi(E2, lead_par.muRu,lead_par.tRu)
        fermi2d  = Fermi(E2,lead_par.muLd,lead_par.tLd) - Fermi(E2, lead_par.muRd,lead_par.tRd)
        
        GfR1     = G_R(sys_par.HS,E1,eta,SR)
        GfR2     = G_R(sys_par.HS,E2,eta,SR)

        GR1_u, GR1_d, GR1_ud, GR1_du = SpinResolvedMatrix(GfR1)
        GR2_u, GR2_d, GR2_ud, GR2_du = SpinResolvedMatrix(GfR2)

        for j in 1:sys_par.N-1
            t1u     = imag.((GR1_u*GLu*(GR1_u)')[j,j+1]*Hu[j+1,j])*fermi1u + imag.((GR1_u*GRu*(GR1_u)')[j,j+1]*Hu[j+1,j])*fermi1u  
            t1d     = imag.((GR1_d*GLd*(GR1_d)')[j,j+1]*Hd[j+1,j])*fermi1d + imag.((GR1_d*GRd*(GR1_d)')[j,j+1]*Hd[j+1,j])*fermi1d  
            t1ud    = imag.((GR1_ud*GLu*(GR1_ud)')[j,j+1]*Hud[j+1,j])*fermi1u + imag.((GR1_ud*GRd*(GR1_ud)')[j,j+1]*Hud[j+1,j])*fermi1d  
            t1du    = imag.((GR1_du*GLu*(GR1_du)')[j,j+1]*Hdu[j+1,j])*fermi1d + imag.((GR1_du*GRu*(GR1_du)')[j,j+1]*Hdu[j+1,j])*fermi1u  
            t2u     = imag.((GR2_u*GLu*(GR2_u)')[j,j+1]*Hu[j+1,j])*fermi2u + imag.((GR2_u*GRu*(GR2_u)')[j,j+1]*Hu[j+1,j])*fermi2u  
            t2d     = imag.((GR2_d*GLd*(GR2_d)')[j,j+1]*Hd[j+1,j])*fermi2d + imag.((GR2_d*GRd*(GR2_d)')[j,j+1]*Hd[j+1,j])*fermi2d  
            t2ud    = imag.((GR2_ud*GLu*(GR2_ud)')[j,j+1]*Hud[j+1,j])*fermi2u + imag.((GR2_ud*GRd*(GR2_ud)')[j,j+1]*Hud[j+1,j])*fermi2d  
            t2du    = imag.((GR2_du*GLu*(GR2_du)')[j,j+1]*Hdu[j+1,j])*fermi2d + imag.((GR2_du*GRu*(GR2_du)')[j,j+1]*Hdu[j+1,j])*fermi2u  
            Inu[j] += dE * (t1u + t2u)
            Ind[j] += dE * (t1d + t2d)
            Inud[j]+= dE * (t1ud + t2ud)
            Indu[j]+= dE * (t1du + t2du)
        end
    end
    return Inu,Ind,Inud,Indu
end

function LB_LocalCurrentV(sys_par,lead_par,Vnu,Vnd,path,savestp)
    mode = "w"
    Inu  = []
    Ind  = []
    Inud = []
    Indu = []
    k    = 0
    for (i,(vu,vd)) in enumerate(zip(Vnu,Vnd))    
        lead_par.muLu = vu/2
        lead_par.muLd = vd/2
        lead_par.muRu = -vu/2
        lead_par.muRd = -vd/2
         
        inu,ind,inud,indu = LB_LocalCurrent(sys_par,lead_par)
        push!(Inu,inu)
        push!(Ind,ind)
        push!(Inud,inud)
        push!(Indu,indu)
         
        if i % savestp == 0
            fil = open(path,mode)
            if k == 0
                write(fil,"#V ")
                for mi in 1:sys_par.N-1
                    write(fil,"Iu$(mi)$(mi+1) ")
                end
                for mi in 1:sys_par.N-1
                    write(fil,"Id$(mi)$(mi+1) ")
                end
                for mi in 1:sys_par.N-1
                    write(fil,"Iud$(mi)$(mi+1) ")
                end
                for mi in 1:sys_par.N-1
                    write(fil,"Idu$(mi)$(mi+1) ")
                end
                
                write(fil,"\n")
            end
            for (v,iu,id,iud,idu) in zip(Vnu[1+k*savestp:(k+1)*savestp],Inu,Ind,Inud,Indu)
                write(fil,"$(v) ")
                for mi in 1:sys_par.N-1
                    write(fil,"$(iu[mi]) ")
                end
                for mi in 1:sys_par.N-1
                    write(fil,"$(id[mi]) ")
                end
                for mi in 1:sys_par.N-1
                    write(fil,"$(iud[mi]) ")
                end
                for mi in 1:sys_par.N-1
                    write(fil,"$(idu[mi]) ")
                end
                write(fil,"\n")
            end
            close(fil)
            mode = "a"
            Inu  = []
            Ind  = []
            Inud = []
            Indu = []
            k += 1
        end
    end
end


end
