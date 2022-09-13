module InputFile

p = expanduser("~/Uni/Master/Masterarbeit/code/src/juliaCode/srcV0.8/")
include("$(p)mainterminal.jl")
include("$(p)hamiltonian.jl")
include("$(p)iosave.jl")
using .MainTerminal
using .Hamiltonian
using .IOsave
using LinearAlgebra

# ---------------- Model -----------------------------

syst = "1d"
leads = false
therm = false

# -------------- Parameters --------------------------

# Simulation Parameters
dt = 0.001
tend = 250
thermstep = Int(1/dt)
numstep = Int(tend/dt)
writeMode = "minimal"
chemTD = "false"
switchMode = "nonSwitched"
hextMode = "tindependent"
aniMode = "none"


# Geometry Parameters
x = 11
y = 1
bc = "(0)"
size = (x,y)
l = Int(x*y)
t = 1.0
PL = []
FL = []
hammeth = "noB"

## Leads
tL = zeros(2)#[1.0,1.0] 
tR = tL# [1.0,1.0] 
V  = 0.0
µ1 = V/2#[V/2 for j in 1:numstep]
µ2 = -V/2
nPoles = 1
beta = 1/0.025

## Central
Jsd = 0.01
Jex = 0.0
B = 0.1
Hext  = [[0.,0.,B] for i in 1:l]
Kani  = [[0.,0.,0.] for i in 1:l]
mag0  = [[1.0,0.0,0.0]]

pos   = [5]

hcen0 = HCentralMulti(l,size,bc,t,Jsd,mag0,pos)
K = 0.0
nfac = 1.0



path = "data"
mkpath(path)
dataname = "$(path)/L$(l)_Jsd$(Jsd)_Jex$(Jex)_B$(B)_K$(K)_V$(V)_dt$(dt)"

params = [("dt",dt),("Ntherm",thermstep),("Nstep",numstep),("chemTD",chemTD),("switchMode",switchMode),("hextMode",hextMode),("aniMode",aniMode),("x",x),("y",y),("bc",bc),("t",t),("tL",tL),("tR",tR),("µ1",µ1),("µ2",µ2),("V",V),("nPoles",nPoles),("beta",beta),("Jsd",Jsd),("Jex",Jex),("mag0",mag0),("Hext",Hext),("Kani",Kani),("nfac",nfac)]

parampath = "$(path)/paramfile.txt"
purpose = "Study of 2d system: TEST. lx = ly = 2, 2 spin in x dir, 2 spin in z dir. Dynamics should be equal to closed system appropriate configuration."
IOsave.SaveParams(params,parampath,purpose)

Simulation(syst,leads,therm,l,Jsd,Jex,mag0,pos,Hext,hcen0,bc,dt,numstep,size,path,hammeth,writeMode,tL,tR,beta,nPoles,µ1,µ2,switchMode,hextMode,aniMode,Kani)

end

