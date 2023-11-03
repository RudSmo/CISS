
import numpy as np
import argparse
import sys
try:
    import scipy.interpolate as interpolate
    import scipy.integrate as integrate
    from scipy import constants as sc
    import matplotlib.pyplot as plt
except:
    pass

from src import mashf90


class CustomArgumentParser(argparse.ArgumentParser):
    """ Allow commented lines in input file """
    def convert_arg_line_to_args(self, line):
        if line.strip().startswith('#'):
            return []
        return line.split()

""" Unit conversions to atomic units """
cmm1 = 4.556335e-6
fs = 41.34136e0
kB = 3.166829e-6
eV = 1./27.2113961


Vconst=0; c=0; nf_site=0
Vlin = 0.; Vbilin = 0.

def read_args():
    parser = CustomArgumentParser(fromfile_prefix_chars=["@","+"])
    conv = parser.add_argument_group("convergence")
    debye = parser.add_argument_group("debye")
    tully = parser.add_argument_group("tully")
    parser.add_argument("-model", type=str, default="spinboson", help="Model system", 
                        choices=["spinboson","biexciton","dba","fmo3","fmo7","fmo8","pyrazine","sf","tully1","tully2","lh2"])
    parser.add_argument("-basis", type=str, default="site", help="Diabatic basis for certain model systems", 
                        choices=["exc","site","adia","dia","rotdia"])
    parser.add_argument("-units", type=str, default="au",help="""Choose unit system. 
                        cmm1: input in cmm1, output in fs. 
                        fs: input in fs, output in fs. """,
                        choices=["au","cmm1","fs"])
    parser.add_argument("-obstyp", type=str, default="pop", help="Observable types", 
                        choices=["pop","all","mannrich","nuc"])
    parser.add_argument("-init", type=int,help="Initial state (in Python indexing)")
    parser.add_argument("-initbasis",type=str,default='dia',choices=["dia","adia","exc"],help='Basis for initial state')
    parser.add_argument("-nucsamp",type=str,default='cl', help="Nuclear sampling (in case 'init' is not None). GS=ground state, WP=wavepacket",
                        choices=['classical','cl','wigner','wig','GS','WP'])
    parser.add_argument("-elsamp",type=str,default='theta',choices=["focused","theta","all","halffoc"],help='Choice of initial distribution')
    parser.add_argument("-cftyp",type=str,default='theta',choices=["cleopatra","wab","theta"],help='Choice of correlation function')
    parser.add_argument("-disorder",type=str,default='none',choices=["none","fmo"],help='Choice of static disorder')
    parser.add_argument("-wo180",action="store_true",help='Remove 180 cm-1 mode in FMO8')
    parser.add_argument("-debug",action="store_true",help='Toggle debug code (plot energy conservation etc)')
    parser.add_argument("-mode",type=str,default='mash', help="Choose which version of SH to run (fssh not fully supported)",
                        choices=['mash','fssh','eql'])
    parser.add_argument("-beta",type=float,default=1.,help="Reciprocal temperature [in a.u.]")
    parser.add_argument("-betamin",type=float,default=1.,help="Reciprocal elevated temperature [in a.u.]")
    parser.add_argument("-T",type=float,default=0,help="Temperature [in kelvin]")
    parser.add_argument("-Tmax",type=float,default=0,help="Maximum temperature for replica exchange")
    parser.add_argument("-useqtb",action="store_true",help="Toggle Quantum Thermal Bath")
    parser.add_argument("-rescalebath",action="store_true",help="Toggle rescaled bath")
    conv.add_argument("-dt",default=41,type=float,help="Time step")
    conv.add_argument("-nt","-TS",type=int,help="Number of time steps")
    conv.add_argument("-ntraj","-traj",type=int,help="Number of trajectories")
    conv.add_argument("-nb",default=1,type=int,help="Number of ring-polymer beads")
    conv.add_argument("-nr",default=1,type=int,help="Number of replicas")
    conv.add_argument("-neql",default=0,type=int,help="Number of equilibriation steps")
    conv.add_argument("-npar","-n",default=1,type=int,help="Number of parallel processors")
    debye.add_argument("-Delta",type=float,default=1.,help="Diabatic coupling")
    debye.add_argument("-epsilon","-eps",type=float,default=0.,help="(Half) energy bias")
    debye.add_argument("-lamda",type=float,default=0,help="System-bath reorganization energy.")
    debye.add_argument("-omegac",type=float,default=0.,help="Cutoff frequency.")
    debye.add_argument("-nf",type=int,default=1,help="Number of bath modes.")
    tully.add_argument("-WPenergy",type=float,default=0.,help="Wavepacket energy")
    tully.add_argument("-gamma",type=float,default=0.,help="Wavepacket width parameter")
    tully.add_argument("-pinit",type=float,default=0.,help="Initial momentum")
    args = parser.parse_args()

    """ ======= Convert input arguments to atomic units ======="""
    if args.units=="cmm1":
        args.epsilon = args.epsilon * cmm1
        args.Delta = args.Delta * cmm1
        args.lamda = args.lamda * cmm1
        args.omegac = args.omegac * cmm1
        args.dt = args.dt * fs
    if args.units=='fs':
        args.dt = args.dt * fs
    if args.T:
        args.beta = 1./(kB*args.T)
    if args.Tmax:
        args.betamin = 1./(kB*args.Tmax)
    
    return args


"""======== Initialize system ========="""
def setup_model(args):
    global Vconst, Vlin, Vbilin, nf_site, c

    model = args.model
    epsilon = args.epsilon
    Delta = args.Delta
    nf = args.nf
    omega = np.zeros(nf,dtype=np.float64)
    if model=='spinboson':
        ns = 2
        Vconst = np.array([[epsilon,Delta],[Delta,-epsilon]])
        Vlin = np.zeros((nf,ns,ns),dtype=np.float64)
        fac = np.sqrt(0.5*args.lamda/nf)
        c = np.zeros(nf)
        for i in range(nf):
            omega[i] = args.omegac*np.tan(0.5*np.pi*(i+0.5)/nf)
            c[i] = fac*omega[i]
            Vlin[i,0,0] = c[i]
            Vlin[i,1,1] = -c[i]
        dim = 0
        print(omega[dim],c[dim])
    elif model=='biexciton':
        ns = 2 
        Vconst = np.array([[0.5*epsilon,Delta],[Delta,-0.5*epsilon]])
        Vlin = np.zeros((nf,ns,ns),dtype=np.float64)
        nf_site = nf//2
        fac = np.sqrt(2*args.lamda/nf_site)
        c = np.zeros(nf_site)
        for i in range(nf_site):
            omega[i] = args.omegac*np.tan(0.5*np.pi*(i+0.5)/nf_site)
            c[i] = fac*omega[i]
        omega[nf_site:] = omega[:nf_site]
        Vlin[:nf_site,0,0] = c
        Vlin[nf_site:,1,1] = c
    elif model=='dba':
        """ State order D B A """
        ns = 3
        Vconst = np.zeros((ns,ns))
        epsB = 20.
        Vconst[1,1] = epsB + args.lamda/4.
        Vconst[2,2] = args.epsilon + args.lamda
        Vconst[0,1] = Vconst[1,0] = Vconst[1,2] = Vconst[2,1] = Delta
        Vlin = np.zeros((nf,ns,ns),dtype=np.float64)
        fac = np.sqrt(2*args.lamda/nf)
        c = np.zeros(nf)
        for i in range(nf):
            omega[i] = args.omegac*np.tan(0.5*np.pi*(i+0.5)/nf)
            c[i] = fac*omega[i]
            Vlin[i,1,1] = 0.5*c[i]
            Vlin[i,2,2] = c[i]
        # print(1./omega)
        # print(Vconst)
    elif model=='fmo3':
        ns = 3
        Vlin = np.zeros((nf,ns,ns),dtype=np.float64)
        omega = np.zeros(nf,dtype=np.float64)
        Vconst = np.diag([12410.,12530,12210])
        Vconst[0,1] = -87.7
        Vconst[0,2] =   5.5
        Vconst[1,2] =  30.8
        Vconst = Vconst + Vconst.T - np.diag(np.diag(Vconst))
        # Vconst -= np.eye(ns)*np.min(np.diag(Vconst))
        Vconst -= np.eye(ns)*np.trace(Vconst)/ns

        Vconst = Vconst*cmm1
        nf_site = nf//ns
        fac = np.sqrt(2*args.lamda/nf_site)
        c = np.zeros(nf_site)
        for i in range(nf_site):
            omega[i] = args.omegac*np.tan(0.5*np.pi*(i+0.5)/nf_site)
            c[i] = fac*omega[i]
        for i in range(ns):
            omega[i*nf_site:(i+1)*nf_site] = omega[:nf_site]
            Vlin[i*nf_site:(i+1)*nf_site,i,i] = -c
    elif model=='fmo7':
        ns = 7
        Vlin = np.zeros((nf,ns,ns),dtype=np.float64)
        omega = np.zeros(nf,dtype=np.float64)
        Vconst = np.diag([12410.,12530,12210,12320,12480,12630,12440])
        Vconst[0,1] = -87.7
        Vconst[0,2] =   5.5
        Vconst[0,3] = - 5.9
        Vconst[0,4] =   6.7
        Vconst[0,5] = -13.7
        Vconst[0,6] = - 9.9
        Vconst[1,2] =  30.8
        Vconst[1,3] =   8.2
        Vconst[1,4] =   0.7
        Vconst[1,5] =  11.8
        Vconst[1,6] =   4.3
        Vconst[2,3] = -53.5
        Vconst[2,4] = - 2.2
        Vconst[2,5] =  -9.6
        Vconst[2,6] =   6.0
        Vconst[3,4] = -70.7
        Vconst[3,5] = -17.0
        Vconst[3,6] = -63.3
        Vconst[4,5] =  81.1
        Vconst[4,6] = - 1.3
        Vconst[5,6] =  39.7
        Vconst = Vconst + Vconst.T - np.diag(np.diag(Vconst))
        Vconst -= np.eye(ns)*np.min(np.diag(Vconst))
        Vconst = Vconst*cmm1
        # import scipy.linalg as linalg
        # rho = linalg.expm(-args.beta*Vconst)
        # print(np.diag(rho)/np.trace(rho))
        # sys.exit()
        nf_site = nf//ns
        fac = np.sqrt(2*args.lamda/nf_site)
        c = np.zeros(nf_site)
        for i in range(nf_site):
            omega[i] = args.omegac*np.tan(0.5*np.pi*(i+0.5)/nf_site)
            c[i] = fac*omega[i]
        for i in range(ns):
            omega[i*nf_site:(i+1)*nf_site] = omega[:nf_site]
            Vlin[i*nf_site:(i+1)*nf_site,i,i] = -c
        
    elif model=="fmo8":
        ns = 8
        Vlin = np.zeros((nf,ns,ns),dtype=np.float64)
        omega = np.zeros(nf,dtype=np.float64)
        Vconst = np.diag([12505.,12425,12195,12375,12600,12515,12465,12700])
        Vconst[0,1] = -94.8
        Vconst[0,2] =   5.5
        Vconst[0,3] = - 5.9
        Vconst[0,4] =   7.1
        Vconst[0,5] = -15.1
        Vconst[0,6] = -12.2
        Vconst[0,7] =  39.5
        Vconst[1,2] =  29.8
        Vconst[1,3] =   7.6
        Vconst[1,4] =   1.6
        Vconst[1,5] =  13.1
        Vconst[1,6] =   5.7
        Vconst[1,7] =   7.9
        Vconst[2,3] = -58.9
        Vconst[2,4] = - 1.2
        Vconst[2,5] =  -9.3
        Vconst[2,6] =   3.4
        Vconst[2,7] =   1.4
        Vconst[3,4] = -64.1
        Vconst[3,5] = -17.4
        Vconst[3,6] = -62.3
        Vconst[3,7] = - 1.6
        Vconst[4,5] =  89.5
        Vconst[4,6] = - 4.6
        Vconst[4,7] =   4.4
        Vconst[5,6] =  35.1
        Vconst[5,7] = - 9.1
        Vconst[6,7] = -11.1
        Vconst = Vconst + Vconst.T - np.diag(np.diag(Vconst))
        Vconst -= np.eye(ns)*np.min(np.diag(Vconst))
        Vconst = Vconst*cmm1

        """ Load datagrabbed spectral density """
        def build_interpolation(x,y):
            spline = interpolate.PchipInterpolator(x,y)
            return lambda x: spline(x) * (x>0)

        Kw_in = np.loadtxt('norm_specden.dat')
        Kw_in[:,1] *= 0.45
        w_cm  = Kw_in[:,0]
        Kw_fs = Kw_in[:,1]
        w_au = w_cm * cmm1
        Kw_au = Kw_fs * fs
        Kw = build_interpolation(w_au,Kw_au)

        Jw = lambda w: Kw(w)*np.pi*w**2

        # reorg = integrate.quad(lambda w: Jw(w)/(np.pi*w), 0, w_au[-1])[0]
        # huangrhys = integrate.quad(lambda w: Jw(w)/(np.pi*w**2), 0, w_au[-1])[0]
        # print('Interpolation moments')
        # print('Huang-Rhys factor: ', huangrhys)
        # print('Reorg. energy: ', reorg/cmm1)

        maxw = 500. * cmm1
        nf_site = nf
        omega_site = maxw * np.arange(1,nf+1)/nf
        Jw_disc = Jw(omega_site)
        dw = omega_site[1]-omega_site[0]
        c = np.sqrt(Jw_disc*omega_site * 2./np.pi*dw)

        if args.wo180:
            w = omega_site
            mask = (w<160*cmm1) | (w>200*cmm1)
            c = np.where(mask,c,0)
            # for w in omega_site:
                # c = np.ma.masked_where(w>160 and w<200,c)
        print('Interpolation moments')
        print('Huang-Rhys factor: ', 0.5*np.sum(c**2/omega_site**3))
        print('Reorg. energy: ', 0.5*np.sum(c**2/omega_site**2)/cmm1)
        # plt.plot(w_au,Jw(w_au),'x')
        # plt.plot(omega_site,Jw_disc/omega_site**2 * dw,'x')
        # plt.plot(omega_site,np.pi/2. * c**2/omega_site**3,'d')
        # plt.title('Jw')
        # plt.figure()
        # plt.plot(omega_site,S_HR,'o')
        # plt.title('SHR')
        # plt.show()
        # sys.exit()

        # """ Load list of Huang-Rhys factors """
        # specden = np.loadtxt('specden.in')
        # omega_site = specden[:,0] * cmm1
        # S_HR = specden[:,1]
        # maxw = 500*cmm1
        # print(omega_site[omega_site<maxw])
        # S_HR = S_HR[omega_site<maxw]
        # omega_site = omega_site[omega_site<maxw]
        # # lamda = specden[:,2] * cmm1
        # c = np.sqrt(2.*omega_site**3*S_HR)
        # Jw = np.pi/2. * c**2/omega_site
        # lamda = np.sum(Jw/omega_site)/np.pi
        # print('S_HR',np.sum(Jw/omega_site**2)/np.pi)
        # print('S_HR',np.sum(S_HR))
        # print('lamda', lamda/cmm1)
        # nf_site = len(omega_site)
        # plt.plot(omega_site,Jw,'x')
        # plt.title('Jw')
        # plt.figure()
        # plt.plot(omega_site,S_HR,'o')
        # plt.title('SHR')
        # plt.show()
        # sys.exit()

        """ Total vibronic part of Hamiltonian """
        nf = ns*nf_site
        omega = np.zeros(nf)
        Vlin = np.zeros((nf,ns,ns))
        for i in range(ns):
            omega[i*nf_site:(i+1)*nf_site] = omega_site
            Vlin[i*nf_site:(i+1)*nf_site,i,i] = c

    elif model=='pyrazine':
        nf = 24
        ns = 2
        import pyrazinef90
        w,ai,bi,ci,aij,bij,cij,bias,rmass=pyrazinef90.pot24()
        omega = w * eV

        delta = bias * eV
        Vconst = np.zeros((ns,ns))
        Vconst[0,0] = -delta
        Vconst[1,1] =  delta
        Vlin =   (np.einsum('i,kl->ikl', ai, np.array([[1,0],[0,0]]))
                + np.einsum('i,kl->ikl', bi, np.array([[0,0],[0,1]]))
                + np.einsum('i,kl->ikl', ci, np.array([[0,1],[1,0]]))) * eV
        Vbilin = (np.einsum('ij,kl->ijkl', aij, np.array(([1,0],[0,0])))
                + np.einsum('ij,kl->ijkl', bij, np.array(([0,0],[0,1])))
                + np.einsum('ij,kl->ijkl', cij, np.array(([0,1],[1,0])))) * eV
    elif model=='sf':
        nf = 28
        ns = 8
        from variables import omega_list, mass_list, linear_coupling, bilin_coupling, stateE_list
        omega = omega_list
        mass = mass_list
        Vconst = stateE_list
        Vlin = linear_coupling
        Vbilin = np.zeros((nf,nf,ns,ns),dtype=np.float64,order='F')
        for i in range(nf):
            Vbilin[i,i] = bilin_coupling[i]
        """ Half bilin to be consistent with what's in the Fortran code """
        Vbilin *= 0.5

        """ TEST: just 2 states """
        ns = 2
        Vconst = Vconst[:ns,:ns]
        Vlin = Vlin[:,:ns,:ns]
        Vbilin = Vbilin[:,:,:ns,:ns]


    elif model=='lh2':
        ns = 24
        """ Ordering as follows: 
        ------B850-----------   ---B800----
        1A 1B 2A 2B ... 8A 8B   1 2 3 ... 8
        """
        """ Load system Hamiltonian """
        Vconst = np.loadtxt('HS.in') * cmm1

        """ Load bath """
        specden = np.loadtxt('specden.in')
        omega_intra = specden[:,0] * cmm1
        S_HR = specden[:,1]/1000.
        lamda = specden[:,2] * cmm1
        c_mono = np.sqrt(2.*omega_intra**3*S_HR)
        Jw = np.pi/2. * c_mono**2/omega_intra
        nf_intra = len(omega_intra)

        """ Ohmic part: J(w) = eta*w*exp(-w/wc)"""
        xi = 0.4
        eta = 2*np.pi*xi # 2*pi/s1^2 = 1
        omegac = 200.*cmm1
        nf_site = args.nf
        w_Ohm = np.zeros(nf_site)
        c_Ohm = np.zeros(nf_site)
        fac = np.sqrt(2*eta*omegac/(nf_site*np.pi))
        for i in range(nf_site):
            w_Ohm[i] = -omegac*np.log((i+0.5)/nf_site)
            c_Ohm[i] = fac*w_Ohm[i]
        Jw_Ohm = np.pi/2.*c_Ohm**2/w_Ohm
        lamda_Ohm = np.sum(Jw_Ohm/w_Ohm)/np.pi
        print(lamda_Ohm/cmm1)

        """ Total bath of one site"""
        nf_site = nf_site + nf_intra
        omega_site = np.concatenate([w_Ohm,omega_intra])
        c = np.concatenate([c_Ohm,c_mono])

        """ Total vibronic part of Hamiltonian """
        nf = ns*nf_site
        omega = np.zeros(nf)
        Vlin = np.zeros((nf,ns,ns))
        for i in range(ns):
            omega[i*nf_site:(i+1)*nf_site] = omega_site
            Vlin[i*nf_site:(i+1)*nf_site,i,i] = c

    """ ======= Initialize mass ======= """
    if 'pyrazine' in model:
        mass = 1./omega
    elif 'sf' in model:
        """ This is handled in variables.py """
        pass
    else:
        mass = np.ones(nf)

    if args.rescalebath:
        """ TEST """
        # x = args.beta*omega/2.
        # f = np.tanh(x)/x
        # for i in range(nf):
        #     Vlin[i] = Vlin[i] * np.sqrt(f[i])
        # print('nf_site',nf_site)
        x = args.beta*omega_site/2.
        f = np.tanh(x)/x
        for i in range(nf_site):
            # print(omega_site[i],f[i])
            c[i] = c[i] * np.sqrt(f[i])
        """ END TEST """

    """ ====== Initialize fortran potential ===== """
    if 'tully' in model:
        nf = 1
        ns = 2
        mass = np.ones(1)*2000. 
        mashf90.init_tully(model[-1],mass)
    elif model in ['pyrazine','sf']:
        mashf90.init_bilinvib(mass,omega,Vconst,Vlin,Vbilin,nf,ns)
    elif model in ['fmo3','fmo7','fmo8','lh2']:
        mashf90.init_frexc(mass,omega,Vconst,c,nf,nf_site,ns)    
    else:
        mashf90.init_linvib(mass,omega,Vconst,Vlin,nf,ns)
    return mass,omega,nf,ns

""" Debugging section to check energy conservation, plot adiabatic populations etc. """
def debug(args,mass,omega,nf,ns):
    beta = args.beta
    dt = args.dt
    nt = args.nt
    """Initialize phase-space variables"""
    q0,p0,qe0,pe0 = sample(args,mass,omega,nf,ns)

    import matplotlib.pyplot as plt
    """ Energy conservation """
    q,p,qe,pe,at,Et,ierr=mashf90.runtrj(q0[:,0], p0[:,0], qe0[:,0], pe0[:,0], dt, nt, nf, ns)
    t = dt*np.arange(nt+1)
    Eref = nf/beta
    plt.plot(t,Et/Eref,'o',alpha=0.5)
    # Ekin = 0.5*np.sum(p**2/mass,1)
    # plt.plot(t,Ekin/Eref,'o',alpha=0.5)
    # plt.figure()
    # plt.plot(t,q)

    """ Adiabatic populations """
    plt.figure()
    qa = np.zeros_like(qe)
    pa = np.zeros_like(pe)
    popa = np.zeros_like(qe)
    for it in range(nt+1):
        qa[it],pa[it] = mashf90.dia2ad(q[it],qe[it],pe[it])
        popa[it] = qa[it]**2 + pa[it]**2
    for n in range(ns):
        label=str(n) 
        plt.plot(t,popa[:,n],'-',color='C%i'%n,alpha=0.5,label=label)
    popat = [popa[it,at[it]-1] for it in range(nt+1)]
    plt.plot(t,popat,'k:',alpha=0.5,label='Active')
    plt.legend()

    if 0:
        dim = 0 
        Tt = np.sum(p**2,1)/2.
        Vt = Et - Tt
        np.savetxt('nuct.out',np.column_stack([t,q[:,dim],p[:,dim]]))
        np.savetxt('Et.out',np.column_stack([t,Et,Vt]))
        np.savetxt('ct.out',np.column_stack([t,qe[:,0],pe[:,0],qe[:,1],pe[:,1]]))
        np.savetxt('cat.out',np.column_stack([t,qa[:,0],pa[:,0],qa[:,1],pa[:,1]]))
    plt.show()

def sample(args,mass,omega,nf,ns):
    """ Sample npar sets of initial phase-space variables.
        qe and pe are real and imaginary parts of c in the diabatic representation """
    npar = args.npar
    beta = args.beta
    q = np.empty((nf,npar),order="F")
    p = np.empty((nf,npar),order="F")
    qe = np.empty((ns,npar),order="F")
    pe = np.empty((ns,npar),order="F")
    if args.disorder=='fmo':
        Vconst_dis = np.copy(Vconst)
        E_sig = np.array([60,100,60,60,120,120,120,100])/(2*np.sqrt(2*np.log(2))) * cmm1
        newE = np.random.normal(np.diag(Vconst), E_sig)
        np.fill_diagonal(Vconst_dis, newE)
        mashf90.init_frexc(mass,omega,Vconst,c,nf,nf_site,ns)    

    for j in range(npar):
        """ Nuclear sampling """
        qmin = 0

        """ TEST: sample around minimum of current state """
        # qmin = -Vlin[:,init,init]/(2*np.diag(Vbilin[:,:,init,init])+mass*omega**2)
        """ END TEST """
        # print(qmin)

        if args.nucsamp=='cl':
            """ Classical """
            q[:,j] = np.random.normal(qmin,1./np.sqrt(beta*mass*omega**2),nf)
            p[:,j] = np.random.normal(0,np.sqrt(mass/beta),nf)
        elif args.nucsamp in ['wigner','wig']:
            """ Wigner """
            qsig = np.sqrt(1./(2*mass*omega*np.tanh(beta*omega/2)))
            psig = np.sqrt(mass*omega/(2*np.tanh(beta*omega/2)))
            q[:,j] = np.random.normal(qmin,qsig,nf)
            p[:,j] = np.random.normal(0,psig,nf)
        elif args.nucsamp=='GS':
            """ Ground state """
            qsig = np.sqrt(1./(2*mass*omega))
            psig = np.sqrt(mass*omega/2.)
            q[:,j] = np.random.normal(qmin,qsig,nf)
            p[:,j] = np.random.normal(0,psig,nf)
        elif args.nucsamp=='WP':
            """ Wave packet """
            if args.gamma>0:
                qsig = 1./np.sqrt(2.*args.gamma)
                psig = np.sqrt(args.gamma/2.)
            else:
                qsig = psig = 0
            qmean = -15.
            if args.WPenergy:
                pmean = np.sqrt(2.*mass*args.WPenergy)
            else:
                pmean = args.pinit
            q[:,j] = np.random.normal(qmean,qsig,nf)
            p[:,j] = np.random.normal(pmean,psig,nf)

        """ Electronic sampling """
        if args.elsamp=='theta':
            """ Theta sampling: sample from region where Theta_init(c) = 1 """
            while 1:
                normals = np.random.randn(2*ns)
                normals = normals/np.linalg.norm(normals)
                qej=normals[:ns]
                pej=normals[ns:]
                pop = qej**2+pej**2
                l=0
                for k in range(ns):
                    if pop[k]>pop[l]: l=k # Find index of state with largest population
                if l==args.init:
                    qe[:,j] = qej; pe[:,j] = pej
                    break
        elif args.elsamp=='focused':
            """ Focussed sampling: sample from points where Phi_n = \delta_{n,init} """
            Hn = np.sum(1./np.arange(1,ns+1))
            A = (ns-1.)/(Hn-1.)
            Pn = np.ones(ns)*(A-1.)/(A*ns)
            Pn[args.init] = 1./ns + (ns-1.)/(A*ns)
            phi = np.random.uniform(0,2*np.pi,ns)
            qe[:,j] = np.cos(phi)*np.sqrt(Pn)
            pe[:,j] = np.sin(phi)*np.sqrt(Pn)
        elif args.elsamp=='halffoc':
            """ half-focussed sampling: sample P_init such that Phi_init = 1 and remaining pops from SU(N-1) sphere """
            Hn = np.sum(1./np.arange(1,ns+1))
            alpha = (ns-1.)/(Hn-1.)
            beta = (1.-alpha)/ns
            Pinit = (1-beta)/alpha
            phi = np.random.uniform(0,2*np.pi)
            qe[args.init,j] = np.cos(phi)*np.sqrt(Pinit)
            pe[args.init,j] = np.sin(phi)*np.sqrt(Pinit)
            normals = np.random.randn(2*ns-2)
            normals = normals/np.linalg.norm(normals) * np.sqrt(1-Pinit)
            qe[np.arange(ns)!=args.init,j] = normals[:ns-1]
            pe[np.arange(ns)!=args.init,j] = normals[ns-1:]
        elif args.elsamp=='all':
            normals = np.random.randn(2*ns)
            normals = normals/np.linalg.norm(normals)
            qe[:,j]=normals[:ns]
            pe[:,j]=normals[ns:]

        """ Perform basis transformation (if requested) """
        if args.initbasis=='adia':
            """ Convert from adia to dia """
            vad,U = mashf90.get_vad(q[:,j],ns)
            qe[:,j] = np.dot(U,qe[:,j])
            pe[:,j] = np.dot(U,pe[:,j])
        elif args.initbasis=='exc':
            """ Convert from exc to dia: transform given by q=0 """
            vad,U = mashf90.get_vad(0.*q[:,j],ns)
            qe[:,j] = np.dot(U,qe[:,j])
            pe[:,j] = np.dot(U,pe[:,j])

    return q,p,qe,pe

def savedata(B,t,ns,args,name):
    tout=t/fs if args.units in ['cmm1','fs'] else t
    if args.obstyp=='pop':
        if 'C' in name:
            shape = (-1,ns**2)
        else:
            shape = (-1,ns)
    elif args.obstyp=='all':
        if 'C' in name:
            shape = (-1,ns**4)
        else:
            shape = (-1,ns**2)
    out = np.column_stack([tout,B.reshape(shape)])
    if name is None: name=args.obstyp
    np.savetxt('%s.out'%name,out)

