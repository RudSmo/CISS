struct HamiltonianParameters
    E1      ::Float64
    E2      ::Float64
    gamma   ::Float64
end

mutable struct PhaseSpaceVariables
    X       ::ComplexF64
    P       ::ComplexF64
    x       ::ComplexF64
    p       ::ComplexF64
    gs      ::Float64
end

struct SimulationParameters
    sphere  ::String
    ntraj   ::Int64
    dt      ::Float64
    Nstp    ::Int64
    Nsave   ::Int64 
    Nflush  ::Int64
    fn      ::String
end

function EOMX(PSVar,HamPar)
    return HamPar.E1 * PSVar.P + HamPar.gamma * PSVar.p
end

function EOMP(PSVar,HamPar)
    return -HamPar.E1 * PSVar.X - HamPar.gamma * PSVar.x 
end

function EOMx(PSVar,HamPar)
    return  HamPar.E2 * PSVar.p + HamPar.gamma * PSVar.P
end

function EOMp(PSVar,HamPar)
    return -HamPar.E2 * PSVar.x - HamPar.gamma * PSVar.X
end


function Populations(PSVar)
    return real.([0.5*(conj(PSVar.X)*PSVar.X + conj(PSVar.P)*PSVar.P - PSVar.gs), 0.5*(conj(PSVar.x)*PSVar.x + conj(PSVar.p)*PSVar.p - PSVar.gs)])
end


function SingleSimulate(PSVar,HamPar,SimPar)
    pops = [Populations(PSVar)]    
    ti   = [0.0]
    for i in 1:SimPar.Nstp
        PSVar.X = PSVar.X + EOMX(PSVar,HamPar)*dt   
        PSVar.P = PSVar.P + EOMP(PSVar,HamPar)*dt
        PSVar.x = PSVar.x + EOMx(PSVar,HamPar)*dt   
        PSVar.p = PSVar.p + EOMp(PSVar,HamPar)*dt
        
        if i % SimPar.Nsave == 0
            push!(pops,Populations(PSVar))
            push!(ti, i*SimPar.dt)
        end
    end
    return ti,pops
end

function NSimulate(HamPar,SimPar,state)
    avg_ind = 2
    PSVar0 = FocusedInitialCondition(SimPar,state)
    ti,Pop_avg = SingleSimulate(PSVar0,HamPar,SimPar)
    for n in 1:SimPar.ntraj-1
        PSVar = FocusedInitialCondition(SimPar,state)
        ti,pops = SingleSimulate(PSVar,HamPar,SimPar)
        Pop_avg .+= pops
        if n % SimPar.Nflush == 0 
            fil = open(SimPar.fn,"w")
            write(fil,"#(1)t (2)P1 (3)P2\n")    
            for j in 1:length(Pop_avg)
                write(fil,"$(ti[j]) $(Pop_avg[j][1]/avg_ind) $(Pop_avg[j][2]/avg_ind)\n")
            end
            close(fil)
        end    
        avg_ind += 1
    end
end

function SphereParams(SimPar)
    """
    Returns Radius of sphere, quaternionic radius, and ZPE parameter
    """
    if SimPar.sphere == "Q"
        rs = 1/2
    elseif SimPar.sphere == "W"
        rs = sqrt(3)/2
    elseif SimPar.sphere == "P"
        rs = 3/2
    else
        error("Sphere not implemented. Try 'Q', 'W', or 'P' instead.") 
    end
    return rs, sqrt(4*rs), 2*rs-1
end

function FocusedInitialCondition(SimPar,state)
    """
    Possible initial states: "S", "B", corresponding to initially occupied system 
    and bath states, respectively.
    """
    rs,Rs,gs = SphereParams(SimPar)
    phiS = rand()*2*pi
    phiB = rand()*2*pi
    if state == "S"
        Js1 = 0
        Js2 = 1
    elseif state == "B"
        Js1 = 1
        Js2 = 0
    else
        error("State not existant. Please use 'S' or 'B' instead.")
    end
    X = sqrt(2*Js1+gs)*cos(phiS)
    P = -sqrt(2*Js1+gs)*sin(phiS)
    x = sqrt(2*Js2+gs)*cos(phiB)
    p = -sqrt(2*Js2+gs)*sin(phiB)     
    return PhaseSpaceVariables(X,P,x,p,gs)
end

E1 = 1.0
E2 = 1.01
gamma = 1.0
sphere = "W"
state = "S"
ntraj = 1001
dt = 0.001
tend = 10
Nstp = Int(tend/dt)
Nsave = 5
Nflush = 9
fn = "data"

if !isdir(fn)
    mkdir(fn)
end


fn = "$(fn)/$(sphere)_E1$(E1)_E2$(E2)_g$(gamma)_ntraj$(ntraj)_dt$(dt).dat"

SimPar = SimulationParameters(sphere,ntraj,dt,Nstp,Nsave,Nflush,fn)
PSVar = FocusedInitialCondition(SimPar,state)
HamPar = HamiltonianParameters(E1,E2,gamma)

NSimulate(HamPar,SimPar,state)

