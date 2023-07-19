import numpy as np
import fnmatch, os
import matplotlib.pyplot as plt


V  = [50,60,65,75,85,105]
Ns = [30,51,51,51,51,31]
ylim = [[0.25,0.5],[0.45,0.85],[0.55,1.05],[0.9,1.4],[1.1,1.65],[1.4,1.7]]
ylim2 = [[-4,4] for i in range(6)]

fig = plt.figure(figsize=(10,14.142))
py = 0.848
px = 0.08
sx = 0.41
sy = 0.1475
dx = 0.08
dy = 0.015

fac_t = 1/(1e3 * 0.002 * 4.135668)
fac_c = 0.002* 1e3 * 0.2433946
size = "2x4"



for iax,v,ns in zip(range(len(V)),V,Ns):

    pxn = px + (sx+dx)
    pyn = py - iax*(sy+dy)


    ax = plt.axes([px,pyn,sx,sy])


    Vfold = f"V{v}/"
    cdatP = np.genfromtxt(Vfold+f"currents_avg_n{ns}_P.dat")
    IP = cdatP[:,1]-cdatP[:,2]
    cdatM = np.genfromtxt(Vfold+f"currents_avg_n{ns}_M.dat")
    IM = cdatM[:,1]-cdatM[:,2]

    ax.plot(cdatP[:,0]*fac_t,IP*fac_c,color="red",linewidth=1)
    ax.plot(cdatM[:,0]*fac_t,IM*fac_c,color="blue",linewidth=1)    

    ax.annotate(r"$V=$"+f"{0.2*0.04*v}"+r"$\,\mathrm{V}$",(0.7,0.75),xycoords="axes fraction",fontsize=13)
    ax.annotate(r"$n_\mathrm{traj}=$"+f"{ns*20}",(0.7,0.6),xycoords="axes fraction",fontsize=13)
    ax.set_ylim(ylim[iax][0],ylim[iax][1])
    ax.set_ylabel(r"$I\,[\mu\mathrm{A}]$")
    if iax != len(V)-1:
        ax.set_xticklabels([])
    else:
        ax.set_xlabel(r"$t\,[\mathrm{ps}]$") 



    ax = plt.axes([pxn,pyn,sx,sy])

    Id = IP-IM
    Is = IP+IM
    SP = 100 * Id / Is


    ax.plot(cdatP[:,0]*fac_t,SP,color="k",linewidth=0.8)
    ax.set_ylabel(r"$M_R\,[\%]$")
    ax.set_ylim(ylim2[iax][0],ylim2[iax][1])
 
    if iax != len(V)-1:
        ax.set_xticklabels([])
    else:
        ax.set_xlabel(r"$t\,[\mathrm{ps}]$")





plt.savefig(f"L{size}_currents_MR_allV_lambda1.pdf")
