import numpy as np
import fnmatch, os
import matplotlib.pyplot as plt


V  = [60,65,75,85]
Ns = [140,131,131,131]
ylim = [[0.55,1.0],[0.7,1.4],[1.1,1.7],[1.4,1.9]]
ylim2 = [[-1.2,1.2],[-1.5,1.5],[-1.5,1.5],[-0.3,2]]

fig = plt.figure(figsize=(10,14.142))
py = 0.77
px = 0.08
sx = 0.41
sy = 0.22
dx = 0.08
dy = 0.015

fac_t = 1/(1e3 * 0.002 * 4.135668)
fac_c = 0.002* 1e3 * 0.2433946
size = "2x2"



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
    ax.annotate(r"$n_\mathrm{traj}=$"+f"{ns*10}",(0.7,0.6),xycoords="axes fraction",fontsize=13)
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





plt.savefig(f"L{size}_currents_MR_allV_longer.pdf")
