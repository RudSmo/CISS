#!/usr/bin/env python3

import numpy as np
import os
import sys

import utils
from runelib import units
from src import mashf90

try:
    import matplotlib.pyplot as plt
except:
    pass

"""=========== Read input parameters=========""" 
args = utils.read_args()

""" ==== Set random seed (comment line if this is not wanted) ==== """
# np.random.seed(42)


""" ======= Save commonly used arguments in their own variables ======"""
model = args.model
beta = args.beta
nf = args.nf
dt = args.dt
nt = args.nt
ntraj = args.ntraj
npar = args.npar
obstyp = args.obstyp

""" Setup paralellization """
if npar>1:
    os.environ["OMP_NUM_THREADS"] = str(npar)

""" ==== Setup time """
t = np.arange(nt+1)*dt

"""======== Initialize potential ========="""
mass, omega, nf, ns = utils.setup_model(args)

if args.mode=='eql':
    """===== Run equilibrium simulation ===="""
    neql = args.neql
    nr = args.nr
    # mf_rpmd = rp.MFrpmd(mass,beta,args.nb,ns)
    # Et,Pt = mf_rpmd.runtrj(dt,nt)
    # nth = int(np.sqrt(neql))
    nth = 0. 
    # nth = 100
    # print(omega[247:255])
    # betamin = beta/2.
    betamin = args.betamin
    # print(beta,betamin)
    mashf90.init_mfrpmd(beta,betamin,args.nb,nr,'p')
    Vt, Tt, Pt = mashf90.run_mfrp(dt,nt+neql,nth,ns,nr)
    # print(Et)
    # t = dt*np.arange(nt+1) / units.fs
    Et = Vt + Tt
    kT = 2*np.mean(Tt[neql:],0)/(args.nb*nf)
    print(1./kT)
    for m in range(nr):
        plt.plot(Et[:,m],alpha=0.3)
        # plt.plot(Vt[:,m],alpha=0.3)
        # plt.plot(Tt[:,m],'--',alpha=0.3)
    plt.plot([neql,neql],[Et.min(),Et.max()],'k--')
    plt.ylabel('Energy')
    plt.tight_layout()
    plt.savefig('Et.png')
    plt.figure()
    t = dt*np.arange(nt+1) / units.fs
    c850 = 'C0'
    c800 = 'C2'
    P850 = np.sum(Pt[neql:,:16],1)
    P800 = np.sum(Pt[neql:,16:],1)
    plt.plot(t,P850,alpha=0.3,color=c850)
    # plt.plot(t,P800,alpha=0.5,color=c800)
    tmp = np.ones_like(t)
    nblock = 10
    lb = nt//nblock
    means = np.mean(Pt[neql:neql+lb*nblock,:].reshape((nblock,-1,ns)),1)
    mean = np.mean(means,0)
    sterr = np.sqrt((np.mean(means**2,0)-mean**2)/(nblock-1))
    # print(np.column_stack([mean,sterr]))
    np.savetxt('Pinf.out',np.column_stack([mean,sterr]))

    means = np.mean(P850[:lb*nblock].reshape((nblock,-1)),1)
    mean = np.mean(means)
    sterr = np.sqrt((np.mean(means**2)-mean**2)/(nblock-1))
    # mean = np.mean(P850)
    # sterr = np.sqrt((np.mean(P850**2)-mean**2)/nt)
    tb = t[:lb*nblock].reshape((nblock,-1))
    tmpb = np.ones_like(tb[0,:])
    for n in range(nblock):
        plt.plot(tb[n,:],tmpb*means[n],'C1--')
    print('P850',mean,sterr)
    plt.fill_between(t,tmp*(mean+2*sterr),tmp*(mean-2*sterr),color=c850,alpha=0.5)
    plt.plot(t,np.ones_like(t)*0.9094,'C3:',label='cl')
    plt.plot(t,np.ones_like(t)*0.90465,'k:',label='quapi')
    plt.title(str(args.nb))
    # means = np.mean(P800[:lb*nblock].reshape((-1,nblock)),0)
    # mean = np.mean(means)
    # sterr = np.sqrt((np.mean(means**2)-mean**2)/(nblock-1))
    # print('P800',mean,sterr)
    # plt.fill_between(t,tmp*(mean+sterr),tmp*(mean-sterr),linestyle='--',color=c800)
    # for n in range(ns):
    #     c = 'C%i'%n
    #     plt.plot(t,Pt[neql:,n],alpha=0.5,color=c)
    #     mean = np.mean(Pt[neql:,n])
    #     sterr = np.sqrt((np.mean(Pt[neql:,n]**2)-mean**2)/nt)
    #     print(n,mean,sterr)
    #     tmp = np.ones_like(t)
    #     plt.fill_between(t,tmp*(mean+sterr),tmp*(mean-sterr),linestyle='--',color=c)
    plt.ylabel('Populations')
    plt.tight_layout()
    plt.legend()
    plt.savefig('Pt.png')
    # plt.figure()
    # plt.hist(St[neql:],bins=100)
    # plt.title('S histogram')
    plt.show()
    sys.exit()

"""===== Initialize MASH Fortran module===="""
mashf90.init_mash(beta, args.useqtb)

""" Debugging section to plot energy conservation, plot adiabatic populations etc. """
if args.debug:
    utils.debug(args,mass,omega,nf,ns)
    sys.exit()

"""Initialize observables"""
if obstyp=='pop':
    Bt = np.zeros((nt+1,ns))
elif obstyp=='all':
    Bt = np.zeros((nt+1,ns,ns),dtype=np.complex128)
if args.elsamp=='all':
    Ct = np.zeros((nt+1,ns,ns))
elif obstyp=='nuc':
    qs = []; ps = []; ws = []
    Tpop = np.zeros((nt+1,ns))
    Rpop = np.zeros((nt+1,ns))

"""======== Loop over trajectories ========"""
ndiscarded = 0
for itraj in range(ntraj//npar):
    q,p,qe,pe = utils.sample(args,mass,omega,nf,ns)

    if obstyp=='pop':
        """ Measure population dynamics """
        if args.elsamp=='all':
            reps = {'site':'d','exc':'e','adia':'a','dia':'d','rotdia':'r'}
            rep = reps[args.basis]
            cftyps = {'cleopatra':'c', 'wab':'w', 'theta':'t'}
            cftyp = cftyps[args.cftyp]
            # args.init only used in case of wab, otherwise set it to something
            if not cftyp=='w' and args.init is None: args.init = 0
            ct, Et, ierr = mashf90.runpar_tcf(q,p,qe,pe,rep,args.init+1,cftyp,dt,nt,nf,ns,npar)
            if args.init is None:
                bt = np.nan
            else:
                bt = ct[:,args.init,:]
        else:
            bt, Et, ierr = mashf90.runpar_poponly(q, p, qe, pe, dt, nt, nf, ns, npar)
    elif obstyp=='all':
        """ Measure dynamics of populations and coherences """
        bt, Et, ierr = mashf90.runpar_all(q, p, qe, pe, dt, nt, nf, ns, npar)    
    elif obstyp=='nuc':
        """ Measure final nuclear distribution (Tully) """
        ierr = np.zeros(npar)
        for j in range(npar):
            qt,pt,qet,pet,at,Et,ierr[j]=mashf90.runtrj(q[:,j], p[:,j], qe[:,j], pe[:,j], dt, nt, nf, ns)
            qs.append(qt[-1,0])
            ps.append(pt[-1,0])

            popt = np.array([mashf90.mash_pops(qt[it],qet[it],pet[it],'d',2) for it in range(nt+1)])
            Tpop += (qt>0)*popt
            Rpop += (qt<0)*popt
            
    """ Check for failed trajectories """
    if sum(ierr)>0:
        ndiscarded += np.sum(ierr>0)
        print('ndiscarded',ndiscarded)
    
    """ Save observables """
    if args.obstyp in ['pop','all']:
        Bt += bt
        if args.elsamp=='all':
            Ct += ct

    """ Store temporary results after each 10 % of the number of trajectories"""
    if ntraj > 10:
        if (itraj+1)%(ntraj//(10*npar)) == 0:
            ctraj = itraj*npar+1 - ndiscarded
            print(ctraj+ndiscarded)
            if args.obstyp in ['pop','all']:
                utils.savedata(Bt/ctraj,t,ns,args,args.obstyp)
                if args.elsamp=='all':
                    utils.savedata(Ct/ctraj,t,ns,args,'C'+args.obstyp)
            np.savetxt('log.out',np.array([ctraj,ntraj]),fmt='%i')

""" Log number of successful trajectories as well as requested number of trajectories """
ctraj = itraj*npar+1 - ndiscarded
np.savetxt('log.out',np.array([ctraj,ntraj]),fmt='%i')
ntraj = ctraj

""" Store final results """
if args.obstyp in ['pop','all']:
    Bt /= ntraj
    utils.savedata(Bt,t,ns,args,args.obstyp)
    if args.elsamp=='all':
        Ct /= ntraj
        utils.savedata(Ct,t,ns,args,'C'+args.obstyp)
if args.obstyp=='nuc':
    bins = 200
    qhist,qbins = np.histogram(qs,bins,density=True)
    phist,pbins = np.histogram(ps,bins,density=True)
    qbins = 0.5*(qbins[:-1]+qbins[1:])
    pbins = 0.5*(pbins[:-1]+pbins[1:])
    np.savetxt('qhist.out',np.column_stack([qbins,qhist]))
    np.savetxt('phist.out',np.column_stack([pbins,phist]))

    Tpop /= ntraj
    Rpop /= ntraj
    np.savetxt('scatt.out',np.column_stack([t/utils.fs,Tpop,Rpop]))
