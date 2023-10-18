p = expanduser("~/Uni/PhD/projects/2022/SOC-Holstein_ComparisonMethods/src/MeanField_NEGF/noSpin/") 
include("$(p)observables.jl")
include("$(p)structures.jl")

using .Observables
using .CustomStructs

src_path=expanduser("~/Uni/PhD/projects/2022/SOC-Holstein_ComparisonMethods/src/MeanField_v3_typed_PVar/")
include("$(src_path)hamiltonian.jl")
include("$(src_path)structures.jl")

t0 = 0.04
eu = -6*t0
ed = -6*t0
E = [eu,ed]
a = 1.0
c = 1.0
tso = t0/40
num_laps = 4
ions = 3
N = Int(num_laps * ions)
tph = zeros(2)
gPh = zeros(2)
tsoph = 0.0
ws = ones(N)*0.002
qs = zeros(N)
model = FranssonHamiltonian
modelStr = "Fransson"
xnIn,pnIn = zeros(N),zeros(N)


Gamma = t0/4

P = -0.5
tl = Gamma
tL = tl * (ones(2) .+ P * [1,-1])
tR = tl * ones(2)

gLu = tL[1]
gLd = tL[2]
gLud = 0.0
gLdu = 0.0

gRu = tR[1]
gRd = tR[2]
gRud = 0.0
gRdu = gRud
HamPar = HamiltonianParams( num_laps,
                            ions,
                            N,
                            model,
                            modelStr,
                            E,
                            a,
                            c,
                            t0,
                            tph,
                            tso,
                            tsoph,
                            gPh,
                            ws,
                            xnIn,
                            pnIn,
                            tL,
                            tR)


HS = model(HamPar)


En = Array(-0.34:0.01:-0.15)
sys_par = SystemParams(N,HS,En)


Temp = 0.025
muL = 0.0
muR = 0.0
tempL  = Temp
tempR  = Temp

lead_par = LeadParamsSF(muL,muR,gLu,gLd,gLud,gLdu,gRu,gRd,gRud,gRdu,tempL,tempR)


Vn = Array(0.0:0.01:1.0)


filname = "fS$(ions)_nL$(num_laps)_tso$(tso)_P$(P)_currents_v3.dat"
savestp = 1
LB_CurrentV(sys_par,lead_par,Vn,filname,savestp)

