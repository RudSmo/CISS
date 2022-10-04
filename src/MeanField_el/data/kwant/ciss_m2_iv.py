#J. Phys. Chem. Lett. 2019, 10, 7126−713
import kwant
import tkwant
import numpy as np
import scipy as sci
import matplotlib.pyplot as plt
import tinyarray
from mpi4py import MPI
import os

#--------------Global variables-----------------
filename='L10Li2so04nomag_v2.dat'
datafile='L10Li2so04nomag_v2'

try:
    os.mkdir(datafile)
except OSError:
    print ("Dir exists")
else:
    print ("Dir OK ")

Ls = 10 #Real system lenght     
Li = 2     #Interfaces lenght 
Lx = Ls + 2*Li   #Total kwant system lenght
Ly = 1

fullC = 10
fullCN = 10

G0 = 5.0/4.0

t = 1.0

Vmax = 6.0*G0

so= G0/5.0

epsL = 0.0 #V-10.0
epsR = 0.0 #V+10.0

muL = epsL
muR = epsR

eps = 0.0


gL = G0/10.0
gR = G0

tL = 3*G0#G0
tR = 3*G0#G0

hL = 0.0

a = 1.0
c = 10.0
e_min = -2.2
e_max =  2.2

ch = 1

Temp=G0/2.0

Tmax = 1.0
NV = 40
NEc=6000
NEt=200


#-------------Global matrices-------------------

sigma_0 = tinyarray.array([[1, 0], [0, 1]])
sigma_u = tinyarray.array([[1, 0], [0, 0]])
sigma_d = tinyarray.array([[0, 0], [0, 1]])
sigma_x = tinyarray.array([[0, 1], [1, 0]])
sigma_y = tinyarray.array([[0, -1j], [1j, 0]])
sigma_z = tinyarray.array([[1, 0], [0, -1]])

#-----------------Lattice and Interfaces--------------
lat = kwant.lattice.square(1,norbs=2)
#There is an extra layer of the leads incorporated into the system
interfaceL = [(lat(Li,0),lat(Li-1,0))]
interfaceR = [(lat(Ls+Li,0),lat(Ls+Li-1,0))]


#-------------------------------------
def vector(i, a, c, fullC):
    phi = i*2*np.pi/fullC
    LL = fullCN 
    return [a*np.cos(phi),a*np.sin(phi),i*c/(LL -1)]

def vectord(i,a,c,fullC,ch):
    [ra,rb,rc] = vector(i, a, c, fullC)
    [r1a,r1b,r1c] = vector(i+ch, a, c, fullC)
    [r2a,r2b,r2c] = vector(i+2*ch, a, c, fullC)
    
    diss1 = np.sqrt((ra-r1a)**2 + (rb-r1b)**2 + (rc-r1c)**2)
    diss2 = np.sqrt((r1a-r2a)**2 + (r1b-r2b)**2 + (r1c-r2c)**2)
    d1 = (np.array([ra,rb,rc]) - np.array([r1a,r1b,r1c]))/diss1
    d2 = (np.array([r1a,r1b,r1c]) - np.array([r2a,r2b,r2c]))/diss2
    print("cross",np.cross(d1,d2))
    return np.cross(d1,d2)
    

def so_sigma_full(i,a,c,fullC,ch):
    [v1,v2,v3] = vectord(i,a, c, fullC,ch)
    return (1j*so*(v1*sigma_x + v2*sigma_y + v3*sigma_z))



def make_system():
#    a,Lx,Ly,t,eps,U,epsL,epsR,tL,tR,gL,gR = params
    system = kwant.Builder()
    for ix in range(Li,Ls+Li):
        iL=ix-Li
        for iy in range(Ly):
            system[lat(ix,iy)] = eps*sigma_0    
            if ix > Li:
                system[lat(ix,iy),lat(ix-1,iy)] = -t*sigma_0 
            if iy > 0:
                system[lat(ix,iy),lat(ix,iy-1)] = -t*sigma_0 
            if ((ix > Li+1) and (ix < Ls+Li)):
                system[lat(ix,iy),lat(ix-2,iy)] = so_sigma_full(iL-2,1.0,c,fullC,ch)
  
        
    # This is here so I can change the voltage drop without redoing the system
    def voltageL(site,V_static):
        return 0.5*V_static*sigma_0 - hL*sigma_z
    
    def voltageR(site,V_static):
        return (-0.5)*V_static*sigma_0
 
    #Few layers of the leads are part of the system
    for ix in range(0,Li):
        for iy in range(Ly):
            system[lat(ix,iy)] = voltageL
            if ix > 0:
                system[lat(ix,iy),lat(ix-1,iy)] = -tL*sigma_0
            if iy > 0:
                system[lat(ix,iy),lat(ix,iy-1)] = -tL*sigma_0
                
    for ix in range(Ls+Li,Lx):
        for iy in range(Ly):
            system[lat(ix,iy)] = voltageR
            if ix > Ls+Li:
                system[lat(ix,iy),lat(ix-1,iy)] = -tR*sigma_0
            if iy > 0:
                system[lat(ix,iy),lat(ix,iy-1)] = -tR*sigma_0
                     
    
    #This accounts for different hopping between leads and the system
    for iy in range(Ly):         
        system[lat(Li,iy),lat(Li-1,iy)] = -gL*sigma_0
        system[lat(Ls+Li,iy),lat(Ls+Li-1,iy)] = -gR*sigma_0
            
    lead_l = kwant.Builder(kwant.TranslationalSymmetry((-1,0)),conservation_law=-sigma_z)
    lead_r = kwant.Builder(kwant.TranslationalSymmetry((1,0)), conservation_law=-sigma_z)        
#    lead_l = kwant.Builder(kwant.TranslationalSymmetry((-1,0)))
#    lead_r = kwant.Builder(kwant.TranslationalSymmetry((1,0)))     
    
    lead_l[(lat(0,iy) for iy in range(Ly))] =   voltageL 
    lead_r[(lat(0,iy) for iy in range(Ly))] =   voltageR
    
    lead_l[lat.neighbors()] = -tL*sigma_0 - 1e-9*sigma_x
    lead_r[lat.neighbors()] = -tR*sigma_0 - 1e-9*sigma_x
    
    system.attach_lead(lead_l)
    system.attach_lead(lead_r)
    return system
#System


dV=Vmax/NV
V=Vmax
system = make_system()
system = system.finalized()
kwant.plot(system)
outfile = open(datafile+"/ciss_dens_"+filename,'w')
outfileu = open(datafile+"/ciss_densu_"+filename,'w')
outfiled = open(datafile+"/ciss_densd_"+filename,'w')

JL = kwant.operator.Current(system, where=interfaceL,sum=True)
JR = kwant.operator.Current(system, where=interfaceR,sum=True)
    
JLu = kwant.operator.Current(system, sigma_u, where=interfaceL,sum=True)
JRu = kwant.operator.Current(system, sigma_u, where=interfaceR,sum=True)
    
JLd = kwant.operator.Current(system, sigma_d, where=interfaceL,sum=True)
JRd = kwant.operator.Current(system, sigma_d, where=interfaceR,sum=True)


D_sys = kwant.operator.Density(system)
    
Du_sys = kwant.operator.Current(system, sigma_u)
Dd_sys = kwant.operator.Current(system, sigma_d)
DSz_sys = kwant.operator.Current(system, sigma_z)

for i1 in range(0,NV+1):
    #def measurement_static(params):
    #kwant.plot(system)
    V=i1*dV
    muL = 0.5*V
    muR = -0.5*V
#    system = make_system()
#    system=system.finalized()
    #Ocupation of leads
    

    
    occupationL = tkwant.manybody.make_occupation(muL,temperature=Temp)
    occupationR = tkwant.manybody.make_occupation(muR,temperature=Temp)
    
    occupation = [occupationL, occupationR]
    # Create the solver
    params = {'V_static': V}
    solver = tkwant.manybody.State(system, Tmax, occupation, params=params)
    solver.refine_intervals(rtol=1E-6, atol=1E-6) 
    #times = range(0.0,Tmax,0.1)

    #for time in times



    solver.evolve(0.5*Tmax)
    resultL = solver.evaluate(JL)
    resultR = solver.evaluate(JR)
    resultLu = solver.evaluate(JLu)
    resultRu = solver.evaluate(JRu)
    resultLd = solver.evaluate(JLd)
    resultRd = solver.evaluate(JRd)
    print("V=",V)
    print("Ls",resultLu,resultLd)
    print("Rs",resultRu,resultRd)
    print("LR",resultL,resultR)
    
#    resD_sys = solver.evaluate(D_sys)
#    resDu_sys = solver.evaluate(Du_sys)
#    resDd_sys = solver.evaluate(D_sys)
#    resDSz_sys = solver.evaluate(DSz_sys)
    
 #   print(V,resD_sys,resDu_sys,resDd_sys,resDSz_sys)
 #   print(V,resultL,resultR,resultLu,resultRu,resultLd,resultRd)
    print(str(V), "\t",str(resultL),"\t",str(resultR),
                  "\t",str(resultLu),"\t",str(resultRu),
                  "\t",str(resultLd),"\t",str(resultRd),file=outfile)
   
    #print(str(V), '   '.join(str(x) for x in resD_sys),file=outfile) 
    #print(str(V), '   '.join(str(x) for x in resDu_sys),file=outfileu)
    #print(str(V), '   '.join(str(x) for x in resDd_sys),file=outfiled)
outfile.close()
outfileu.close()
outfiled.close()

#ennp = np.array([])
#daLR = np.array([])
#daRL = np.array([])
#daLRuu = np.array([])
#daLRdd = np.array([])
#daLRud = np.array([])
#daLRdu = np.array([])
##J_0 = kwant.operator.Current(system)
#
#params = {'V_static': V}
#
#E_min =-2.0*(max(tL,tR))-0.5*V+0.00001
#E_max =2.0*(max(tL,tR))+0.5*V-0.00001
#de = (e_min - E_min)/(NEt)
#ennp = np.array([])
#for ie in range(NEt):
#    ennp = np.append(ennp, E_min + de*ie)
#
#de = (e_max - e_min)/(NEc)    
#for ie in range(NEc):
#    energy = e_min + ie*de
#    ennp = np.append(ennp,energy)
# 
#de = (E_max - e_max)/(NEt)    
#for ie in range(NEt+1):
#    energy = e_max + ie*de
#    ennp = np.append(ennp,energy)
#    
#
#rho_sz = kwant.operator.Density(system,sigma_z)
#rho_all = kwant.operator.Density(system)
#
#rho_u = kwant.operator.Density(system,sigma_u)
#rho_d = kwant.operator.Density(system,sigma_d)
#
#rho_Sz = []
#rho_All = []
#rho_U = []
#rho_D = []
#for energy in ennp:
##    print(energy)
##    energy = e_min + de*ie
#    smatrix = kwant.smatrix(system,energy,params=params)
##    ennp = np.append(ennp,energy)
#    daLR=np.append(daLR,smatrix.transmission(1,0))
#    daRL=np.append(daRL,smatrix.transmission(0,1))
#    daLRuu=np.append(daLRuu,smatrix.transmission((1,0),(0,0)))
#    daLRdd=np.append(daLRdd,smatrix.transmission((1,1),(0,1)))
#    daLRud=np.append(daLRud,smatrix.transmission((1,0),(0,1)))
#    daLRdu=np.append(daLRdu,smatrix.transmission((1,1),(0,0)))
#    rho_szb = rho_sz.bind(params=params)
#    rho_allb = rho_all.bind(params=params)
#    rho_ub = rho_u.bind(params=params)
#    rho_db = rho_d.bind(params=params)
#    wf = kwant.wave_function(system, energy, params=params)
#    all_w = np.vstack((wf(0), wf(1)))
#    rho_Sz.append(sum(rho_szb(p) for p in all_w))
#    rho_All.append(sum(rho_allb(p) for p in all_w))
#    rho_U.append(sum(rho_ub(p) for p in all_w))
#    rho_D.append(sum(rho_db(p) for p in all_w))
#
#
#rhoSz = np.array(rho_Sz)
#rhoAll = np.array(rho_All)
#rhoU = np.array(rho_U)
#rhoD = np.array(rho_D)
#
#plt.figure()
#plt.plot(ennp,daLR,'b.')
#plt.plot(ennp,daRL,'r-')
#plt.plot(ennp,daLRuu,'y.')
#plt.plot(ennp,daLRdd,'g-')
#plt.plot(ennp,daLRud,'k.')
#plt.plot(ennp,daLRdu,'p-')
#
#plt.plot(ennp,daLRuu-daLRdd,'r-')
#
##plt.plot(ennp,data3np-data4np,'b.')
#
#plt.xlabel("energy [t]")
#plt.ylabel("conductance [e^2/h]")
#plt.xlim((e_min, e_max))
#plt.show()
#
#plt.plot(ennp,daLRuu-daLRdd,'r-')
#plt.xlim((e_min, e_max))
#plt.show()
#
#rhoSzE = (sum(rhoSz[:,i] for i in range(Li,Ls+Li)))
#rhoAllE = (sum(rhoAll[:,i] for i in range(Li,Ls+Li)))
#rhoUE = (sum(rhoU[:,i] for i in range(Li,Ls+Li)))
#rhoDE = (sum(rhoD[:,i] for i in range(Li,Ls+Li)))
#
#
#plt.plot(ennp,rhoSzE,'r.')
#plt.plot(ennp,rhoAll[:,0],'r-')
#plt.plot(ennp,rhoAll[:,6],'b-')
#
#plt.plot(ennp,rhoAll[:,6],'r-')
#plt.xlim((e_min, e_max))
#plt.show()
#
#plt.plot(ennp,rhoAllE,'b.')
#plt.xlim((e_min, e_max))
#plt.show()
#
#plt.plot(ennp,rhoUE,'r.')
#plt.xlim((e_min, e_max))
#plt.show()
#
#plt.plot(ennp,rhoDE,'b.')
#plt.xlim((e_min, e_max))
#plt.show()
#
#plt.plot(ennp,rhoUE-rhoDE,'b-')
#plt.xlim((e_min, e_max))
#plt.show()
#
#rhoSzL=np.array([])
#rhoAllL=np.array([])
#IL = np.array([])
#
#
#  
#for iN in range(Lx):
##    rhoSzL=np.append(rhoSzL,np.trapz(rhoSz[:,iN],ennp))
##    rhoAllL=np.append(rhoAllL,np.trapz(rhoAll[:,iN],ennp))
#    rhoSzL=np.append(rhoSzL,sci.integrate.simps(rhoSz[:,iN],ennp))
#    rhoAllL=np.append(rhoAllL,sci.integrate.simps(rhoAll[:,iN],ennp))
#    IL = np.append(IL,iN)
#
#
#plt.plot(IL,rhoSzL,'b-')
#plt.show()
#plt.plot(IL,rhoAllL/(2*np.pi),'b*')
#plt.show()
#
#outfile = open(datafile+"/ciss_en_trans_"+filename,'w')
#print("#en,TLR,TRL,Tuu,Tdd,Tud,Tdu",file=outfile)
#for i1 in range(len(ennp)):
#       print(str(ennp[i1]), "\t",str(daLR[i1]),"\t",
#            str(daRL[i1]),"\t",str(daLRuu[i1]),"\t",
#            str(daLRdd[i1]),"\t",str(daLRud[i1]),"\t",
#            str(daLRdu[i1]),file=outfile) 
#outfile.close()
#
#outfile = open(datafile+"/ciss_en_dens_"+filename,'w')
#print("#en,rho,rho_Sz,rho_u,rho_d",file=outfile)
#for i1 in range(len(ennp)):
#       print(str(ennp[i1]), "\t",str(rhoSzE[i1]),"\t",
#            str(rhoAllE[i1]),"\t",str(rhoUE[i1]),"\t",
#            str(rhoDE[i1]),file=outfile) 
#outfile.close()
#
#outfile = open(datafile+"/ciss_en_densSzL_"+filename,'w')
#print("#en,rho_Sz",file=outfile)
#for i1 in range(len(ennp)):
#       print(str(ennp[i1]), '   '.join(str(x) for x in rhoSz[i1,:]),file=outfile) 
#outfile.close()
#
#outfile = open(datafile+"/ciss_en_densL_"+filename,'w')
#print("#en,rho_Sz",file=outfile)
#for i1 in range(len(ennp)):
#       print(str(ennp[i1]), '   '.join(str(x) for x in rhoAll[i1,:]),file=outfile) 
#outfile.close() 
#
#outfile = open(datafile+"/ciss_en_densUL_"+filename,'w')
#print("#en,rho_Sz",file=outfile)
#for i1 in range(len(ennp)):
#       print(str(ennp[i1]), '   '.join(str(x) for x in rhoU[i1,:]),file=outfile) 
#outfile.close()
#
#outfile = open(datafile+"/ciss_en_densDL_"+filename,'w')
#print("#en,rho_Sz",file=outfile)
#for i1 in range(len(ennp)):
#       print(str(ennp[i1]), '   '.join(str(x) for x in rhoD[i1,:]),file=outfile) 
#outfile.close()
#
#
##energy_subset = np.linspace(-2, 2)
##spectrum = kwant.kpm.SpectralDensity(system,params=params)
##densities = spectrum(energy_subset)
#
##plt.plot(energy_subset,densities)
##plt.plot(IL,rho0,'b.')rhoAllL/(2*np.pi)