clear
close all
clc

# Hamiltonian Parameters
num_laps = 1;
num_ions = 3;
N = num_laps * num_ions;
ws  = 0.002*ones(N);
T0  = 0.04/ws(1) ; # in eV
t0  = -T0;
eu  = -6*T0;
ed  = -6*T0;
E   = [eu,ed];
gPh = zeros(2);
a   = 1.0;
c   = 1.0;
tso = T0/40;
tPh = 0.1*T0*ones(2);
tsoPh = 0.1*tso;



#xnIn,pnIn = [0.7691068624482689,-3.61069243320506,-5.245657309206734],[-7.401418656594664,-4.3318963287495675,1.5952234947490966]
xnIn = zeros(N);
pnIn = zeros(N);

## Leads
P = 0.5;
tl = T0/4;
tL = tl * (ones(2) .+ P * [1,-1])  ;
tR = tl * ones(2);

ham = FranssonHamiltonian(N,E,t0,tso,tPh,tsoPh,a,c,num_laps,ws,xnIn)





