module MappingNLSPhonons

using LinearAlgebra
using ..Structures: PhaseSpaceVariables
export NSimulate,FocusedInitialCondition

function EOMX(PSVar,HamPar,SimPar,ip)
    if ip % 2 == 0
        n = Int(ip/2)
    else
        n = Int((ip+1)/2)
    end
    return sum([2 * HamPar.Hamiltonian[ip,j] * PSVar.P[j] for j in 1:SimPar.Nlvl]) + sqrt(2) * HamPar.gamma[n]*PSVar.P[ip] * PSVar.x[n]
end

function EOMP(PSVar,HamPar,SimPar,ip)
    if ip % 2 == 0
        n = Int(ip/2)
    else
        n = Int((ip+1)/2)
    end
    return -sum([2 * HamPar.Hamiltonian[ip,j] * PSVar.X[j] for j in 1:SimPar.Nlvl]) - sqrt(2) * HamPar.gamma[n]*PSVar.X[ip] * PSVar.x[n]
end

function EOMx(PSVar,HamPar,n)
    return HamPar.omega[n] * PSVar.p[n]    
end

function EOMp(PSVar,HamPar,ip)
    if ip % 2 == 0
        n = Int(ip/2)
        k = -1
    else
        n = Int((ip+1)/2)
        k = 1
    end
    return -HamPar.omega[n] * PSVar.x[n] - HamPar.gamma[n]/sqrt(2) * ((PSVar.X[ip]^2 + PSVar.P[ip]^2 - PSVar.gs) + (PSVar.X[ip+k]^2 + PSVar.P[ip+k]^2 - PSVar.gs))
end


function Populations(PSVar,SimPar)
    return real.([0.5*(conj(PSVar.X[i])*PSVar.X[i] + conj(PSVar.P[i])*PSVar.P[i] - PSVar.gs) for i in 1:SimPar.Nlvl])
end

function SingleSimulate(PSVar,HamPar,SimPar)
    pops = [Populations(PSVar,SimPar)]    
    xl   = [PSVar.x]
    pl   = [PSVar.p]
    ti   = [0.0]
    for i in 1:SimPar.Nstp
        for j in 1:SimPar.Nlvl
            PSVar.X[j] = PSVar.X[j] + EOMX(PSVar,HamPar,SimPar,j)*SimPar.dt   
            PSVar.P[j] = PSVar.P[j] + EOMP(PSVar,HamPar,SimPar,j)*SimPar.dt
            if j % 2 != 0
                n = Int((j+1)/2)
                PSVar.x[n] = PSVar.x[n] + EOMx(PSVar,HamPar,n)*SimPar.dt   
                PSVar.p[n] = PSVar.p[n] + EOMp(PSVar,HamPar,j)*SimPar.dt
            end
        end
        if i % SimPar.Nsave == 0
            push!(pops,Populations(PSVar,SimPar))
            push!(xl,real(PSVar.x))
            push!(pl,real(PSVar.p))
            push!(ti, i*SimPar.dt)
        end
    end
    return ti,pops,xl,pl
end

function NSimulate(HamPar,SimPar,state)
    avg_ind = 2
    PSVar0 = FocusedInitialCondition(SimPar,HamPar,state)
    ti,Pop_avg,xl_avg,pl_avg = SingleSimulate(PSVar0,HamPar,SimPar)
    for n in 1:SimPar.ntraj-1
        PSVar = FocusedInitialCondition(SimPar,HamPar,state)
        ti,pops,xl,pl = SingleSimulate(PSVar,HamPar,SimPar)
        Pop_avg .+= pops
        xl_avg  .+= xl
        pl_avg  .+= pl
        if n % SimPar.Nflush == 0 
            #rm(SimPar.fn)
            filP = open(SimPar.fnP,"w")
            filE = open(SimPar.fnE,"w")
            write(filE,"#(1)t ")    
            write(filP,"#(1)t ")
            for j in 1:SimPar.Nlvl
                write(filE,"#($(j+1))P$(j) ")
                if j % 2 != 0
                    n = Int((j+1)/2)
                    write(filP,"#($(2*n-1))x$(n) ($(2*n))p$(n) ")
                end
            end
            write(filE,"\n")
            write(filP,"\n")
            for j in 1:length(Pop_avg)
                write(filE,"$(ti[j]) ")
                write(filP,"$(ti[j]) ")
                for k in 1:SimPar.Nlvl
                    write(filE,"$(Pop_avg[j][k]/avg_ind) ")
                    if k % 2 != 0
                        n = Int((k+1)/2)
                        write(filP,"$(real(xl_avg[j][n])/avg_ind) $(real(pl_avg[j][n])/avg_ind) ")
                    end
                end
                write(filE,"\n")
                write(filP,"\n")
            end
            close(filE)
            close(filP)
        end    
        avg_ind += 1
    end
end

function SphereParams(SimPar)
    """
    Returns Radius of sphere, quaternionic radius, and ZPE parameter
    """
    if SimPar.sphere == "Q"
        Rs = 2
    elseif SimPar.sphere == "W"
        Rs = 2*sqrt(SimPar.Nlvl+1)
    elseif SimPar.sphere == "P"
        Rs = 2*(SimPar.Nlvl+1)
    else
        error("Sphere not implemented. Try 'Q', 'W', or 'P' instead.") 
    end
    return Rs, (Rs^2-2)/SimPar.Nlvl
end

function FocusedInitialCondition(SimPar,HamPar,state)
    """
    Possible initial states: "S", "B", corresponding to initially occupied system 
    and bath states, respectively.
    """
    Rs,gs = SphereParams(SimPar)
    phis = rand(SimPar.Nlvl)*2pi
    X = [sqrt(2*state[j]+gs)*cos(phis[j]) for j in 1:SimPar.Nlvl]
    P = [-sqrt(2*state[j]+gs)*sin(phis[j]) for j in 1:SimPar.Nlvl]
    x = []
    p = []
    for n in 1:Int(SimPar.Nlvl/2)    
        w  = HamPar.omega[n]
        aj = tanh(0.5*HamPar.beta*w)
        A  = randn(2)*aj/pi
        push!(x,A[1]/sqrt(w*aj))
        push!(p,A[2]*sqrt(w/aj))
    end
    return PhaseSpaceVariables(X,P,x,p,gs)
end

end