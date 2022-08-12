struct HamiltonianParameters
    """
    ES = System energies
    EB = Bath energies
    gamma = Couplings

    Ordering (Example): Three-level system
    E = [ES1,ES2,EB1]
    gamma = [[gS1S2,gS1B1],[gS2S1,gS2B1],[gB1S1,gB1S2]]
    """
    Hamiltonian ::Array{Float64}
end

mutable struct PhaseSpaceVariables
    """
    Ordering (Example): Three-level system
    X = [XS1,XS2,XB1]
    P = [PS1,PS2,PB1]
    """
    X       ::Array{ComplexF64}
    P       ::Array{ComplexF64}
    gs      ::Float64
end

struct SimulationParameters
    Nlvl    ::Int64
    sphere  ::String
    ntraj   ::Int64
    dt      ::Float64
    Nstp    ::Int64
    Nsave   ::Int64 
    Nflush  ::Int64
    fn      ::String
end

function EOMX(PSVar,HamPar,ip)
    return sum([HamPar.Hamiltonian[ip,j] * PSVar.P[j] for j in 1:SimPar.Nlvl])
end

function EOMP(PSVar,HamPar,ip)
    return -sum([HamPar.Hamiltonian[ip,j] * PSVar.X[j] for j in 1:SimPar.Nlvl])
end


function Populations(PSVar)
    return real.([0.5*(conj(PSVar.X[i])*PSVar.X[i] + conj(PSVar.P[i])*PSVar.P[i] - PSVar.gs) for i in 1:SimPar.Nlvl])
end


function SingleSimulate(PSVar,HamPar,SimPar)
    pops = [Populations(PSVar)]    
    ti   = [0.0]
    for i in 1:SimPar.Nstp
        for j in 1:SimPar.Nlvl
            PSVar.X[j] = PSVar.X[j] + EOMX(PSVar,HamPar,j)*dt   
            PSVar.P[j] = PSVar.P[j] + EOMP(PSVar,HamPar,j)*dt
        end
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
            #rm(SimPar.fn)
            fil = open(SimPar.fn,"w")
            write(fil,"#(1)t \n")    
            for j in 1:SimPar.Nlvl
                write(fil,"#($(j+1))P$(j) ")
            end
            write(fil,"\n")
            for j in 1:length(Pop_avg)
                write(fil,"$(ti[j]) ")
                for k in 1:SimPar.Nlvl
                    write(fil,"$(Pop_avg[j][k]/avg_ind) ")
                end
                write(fil,"\n")
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

function FocusedInitialCondition(SimPar,state)
    """
    Possible initial states: "S", "B", corresponding to initially occupied system 
    and bath states, respectively.
    """
    Rs,gs = SphereParams(SimPar)
    phis = rand(SimPar.Nlvl)*2pi
    X = [sqrt(2*state[j]+gs)*cos(phis[j]) for j in 1:SimPar.Nlvl]
    P = [-sqrt(2*state[j]+gs)*sin(phis[j]) for j in 1:SimPar.Nlvl]
    return PhaseSpaceVariables(X,P,gs)
end

Nlvl = 3
ES1 = 1.0
ES2 = -1.0
EB1 = 0.0
E = [ES1,ES2,EB1]
gsb = 1.0
Hamiltonian = [ES1 0.0 gsb;
               0.0 ES2 gsb;
               gsb gsb EB1]


sphere = "W"
state = [1,0,0]

ntraj = 1001
dt = 0.001
tend = 10
Nstp = Int(tend/dt)
Nsave = 5
Nflush = 9
fn = "data_3LS"

if !isdir(fn)
    mkdir(fn)
end


fn = "$(fn)/$(sphere)_S1$(ES1)_S2$(ES2)_EB$(EB1)_gsb$(gsb)_ntraj$(ntraj)_dt$(dt).dat"

SimPar = SimulationParameters(Nlvl,sphere,ntraj,dt,Nstp,Nsave,Nflush,fn)
PSVar = FocusedInitialCondition(SimPar,state)
HamPar = HamiltonianParameters(Hamiltonian)

NSimulate(HamPar,SimPar,state)

