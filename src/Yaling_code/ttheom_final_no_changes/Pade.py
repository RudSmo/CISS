import numpy as np

class Pade_decomposition():

    def __init__(self,number_Pade_poles,Spec):
        self.number_Pade_poles=number_Pade_poles
        self.Spec=Spec
        if (self.number_Pade_poles > 0):
            if 'Fermi' in self.Spec or 'fermi' in self.Spec:
                self.Generate_Fermi_Pade_Parameters()
            elif 'Bose' in self.Spec or 'bose' in self.Spec:
                self.Generate_Bose_Pade_Parameters()


    def Generate_Fermi_Pade_Parameters(self):
        """
        Function generating the Pade poles and their respective weights. The routine used here is the [N-1,N] PSD scheme proposed in JCP,134,244106(2011).
        """
        def Kronecker_Delta(m,n):
            if m==n:
                return 1.
            else:
                return 0.
        def b_Fermi(m):
            return 2.*m-1.
    
        M=2*self.number_Pade_poles
        Delta       =np.zeros((M,M),float)
        Delta_tilde =np.zeros((M-1,M-1),float) 
    
        for m in range(M):
            for n in range(M):
                Delta[m,n] =(Kronecker_Delta(m+1,n+1+1)+Kronecker_Delta(m+1,n+1-1))/np.sqrt( b_Fermi(m+1)*b_Fermi(n+1))
    
        for m in range(M-1):
            for n in range(M-1):
                Delta_tilde[m,n] =(Kronecker_Delta(m+1,n+1+1)+Kronecker_Delta(m+1,n+1-1))/np.sqrt( b_Fermi(m+2)*b_Fermi(n+2))
    
        roots_Delta         =np.linalg.eigvalsh(Delta)
        roots_Delta_tilde   =np.linalg.eigvalsh(Delta_tilde)
    
        self.Xis_Fermi=[]
        for root in roots_Delta:
            if root>0:
                self.Xis_Fermi.append(2./abs(root))
        self.Xis_Fermi.sort()
        self.Chis_Fermi=[]
        for root in roots_Delta_tilde:
            if abs(root)>=1.e-13:
                if root>0:
                    self.Chis_Fermi.append(2./abs(root))
        self.Chis_Fermi.sort()
    
        self.kappas_Fermi=[]
        for j in range(self.number_Pade_poles):
            Quotient=1.
            for k in range(self.number_Pade_poles):
                if k!=(self.number_Pade_poles-1) and k!=j:
                    Quotient*=(self.Chis_Fermi[k]**2-self.Xis_Fermi[j]**2)/(self.Xis_Fermi[k]**2-self.Xis_Fermi[j]**2)
                elif k==(self.number_Pade_poles-1) and k!=j:
                    Quotient*=1./(self.Xis_Fermi[k]**2-self.Xis_Fermi[j]**2)
                elif k!=(self.number_Pade_poles-1) and k==j:
                    Quotient*=(self.Chis_Fermi[k]**2-self.Xis_Fermi[j]**2)
            Quotient*=0.5*self.number_Pade_poles*b_Fermi(self.number_Pade_poles+1)
            self.kappas_Fermi.append(Quotient)
    
    def Generate_Bose_Pade_Parameters(self):
        """
        Function generating the Pade poles and their respective weights. The routine used here is the [N-1,N] PSD scheme proposed in JCP,134,244106(2011).
        """
        def Kronecker_Delta(m,n):
            if m==n:
                return 1.
            else:
                return 0.
        
        def b_Bose(m):
            return 2.*m+1.
    
        M=2*self.number_Pade_poles
        Delta       =np.zeros((M,M),float)
        Delta_tilde =np.zeros((M-1,M-1),float) 
        for m in range(M):
            for n in range(M):
                Delta[m,n] =(Kronecker_Delta(m+1,n+1+1)+Kronecker_Delta(m+1,n+1-1))/np.sqrt( b_Bose(m+1)*b_Bose(n+1))
    
        for m in range(M-1):
            for n in range(M-1):
                Delta_tilde[m,n] =(Kronecker_Delta(m+1,n+1+1)+Kronecker_Delta(m+1,n+1-1))/np.sqrt( b_Bose(m+2)*b_Bose(n+2))
    
        roots_Delta         =np.linalg.eigvalsh(Delta)
        roots_Delta_tilde   =np.linalg.eigvalsh(Delta_tilde)
    
        self.Xis_Bose=[]
        for root in roots_Delta:
            if root >0:
                self.Xis_Bose.append(2./root)
        self.Xis_Bose.sort()
        self.Chis_Bose=[]
        for root in roots_Delta_tilde:
            if abs(root)>=1.e-13:
                if root >0:
                    self.Chis_Bose.append(2./abs(root))
        self.Chis_Bose.sort()
    
        self.kappas_Bose=[]
        for j in range(self.number_Pade_poles):
            Quotient=1.
            for k in range(self.number_Pade_poles):
                if k!=(self.number_Pade_poles-1) and k!=j:
                    Quotient*=(self.Chis_Bose[k]**2-self.Xis_Bose[j]**2)/(self.Xis_Bose[k]**2-self.Xis_Bose[j]**2)
                elif k==(self.number_Pade_poles-1) and k!=j:
                    Quotient*=1./(self.Xis_Bose[k]**2-self.Xis_Bose[j]**2)
                elif k!=(self.number_Pade_poles-1) and k==j:
                    Quotient*=(self.Chis_Bose[k]**2-self.Xis_Bose[j]**2)
            Quotient*=0.5*self.number_Pade_poles*b_Bose(self.number_Pade_poles+1)
            self.kappas_Bose.append(Quotient)
    
    
    def Fermi_Pade_Approx(self, energy=0.,mu=0.,temperature=20.):
        """
        Approximation of the Fermi-Dirac distribution function using the Pade poles 'etas' and 'Xis'.
        Input: energy, chemical potential mu, and temperatue are in units of eV. Pade poles, etas and Xis, are dimensionless.
        Output: approximate Fermi-Dirac distribution function.
        """
        f=0.5
        x=(energy-mu)/temperature
        for j in range(self.number_Pade_poles):
            f-=2.*self.kappas_Fermi[j]*x/(x**2+self.Xis_Fermi[j]**2)
        return f
    
    def Fermi_V(self,energy=0.,mu=0.,temperature=20.):
        """
        Input: energy, the chemical potential mu, and the temperature are in units of eV.
        Output: Fermi-Dirac Distribution Function.
        """
        try:
            return 1./( 1.+np.exp( (energy-mu)/temperature ) )
        except:
            return 0.
    
    def Bose_Pade_Approx(self,energy=0.,temperature=10.):
        """
        Approximation of the Bose-Einstein distribution function using the Pade poles 'etas' and 'Xis'.
        Input: energy, chemical potential mu, and temperatue are in units of eV, eV, and K, respectively. Pade poles, etas and Xis, are dimensionless.
        Output: approximate Bose-Einstein distribution function.
        """
        x=energy/temperature
        f=0.5+1./x
        for j in range(self.number_Pade_poles):
            f+=2.*self.kappas_Bose[j]*x/(x**2+self.Xis_Bose[j]**2)
        return f

    def Bose_V(self,energy=0.1,temperature=20.):
        """
        Input: energy and the temperature are in units of eV.
        Output: Bose-Einstein Distribution Function.
        """
        try:
            return 1./( 1.-np.exp( -energy/temperature ) )
        except:
            return 0.

    def get_Pade_parameters(self):
        if 'Fermi' in self.Spec or 'fermi' in self.Spec:
            return self.kappas_Fermi,self.Xis_Fermi
        elif 'Bose' in self.Spec or 'bose' in self.Spec:
            return self.kappas_Bose,self.Xis_Bose

    def get_coefficients_exponents(self, temperature, broadenings, number_bath, number_spin = 1, chem_potentials = 0., WBL=False):
        if 'Bose' in self.Spec or 'bose' in self.Spec:
            coefficients=np.zeros((number_bath,(self.number_Pade_poles+1)),complex)
            exponents=np.zeros((number_bath,(self.number_Pade_poles+1)),complex)
            Deltas=np.zeros((number_bath),float)
            for i_bath in range(number_bath):
                Deltas[i_bath] += 2./broadenings[i_bath] * temperature[i_bath] - 1./np.tan(broadenings[i_bath]/(2.0*temperature[i_bath]))
                for i_pole in range(self.number_Pade_poles+1):
                    if i_pole==0:
                        coefficients[i_bath,i_pole]=broadenings[i_bath]*(1./np.tan(broadenings[i_bath]/(2.0*temperature[i_bath]))-1.j)
                        #coefficients[i_bath,i_pole]=broadenings[i_bath]*(2.0*temperature[i_bath]/broadenings[i_bath]-1.j)
                        exponents[i_bath,i_pole]=broadenings[i_bath]
                    else:
                       #coefficients[i_bath,i_pole]=8.*i_pole*np.pi*broadenings[i_bath]/((2.*i_pole*np.pi)**2-(broadenings[i_bath]/temperature[i_bath])**2)      
                       #exponents[i_bath,i_pole]=2.*i_pole*np.pi*temperature[i_bath]
                        coefficients[i_bath,i_pole]=4.*broadenings[i_bath]*self.Xis_Bose[i_pole-1]*self.kappas_Bose[i_pole-1]/((self.Xis_Bose[i_pole-1])**2-(broadenings[i_bath]/temperature[i_bath])**2)      
                        exponents[i_bath,i_pole]=self.Xis_Bose[i_pole-1]*temperature[i_bath]
                        Deltas[i_bath] -= np.real( (coefficients[i_bath,i_pole] / exponents[i_bath,i_pole]))
            return Deltas, coefficients, exponents

        elif 'Fermi' in self.Spec or 'fermi' in self.Spec:
            if WBL == True:
                coefficients=np.zeros((number_bath, number_spin, self.number_Pade_poles),complex)
                for i_bath in range(number_bath):
                    for i_spin in range(number_spin):
                        for i_pole in range(self.number_Pade_poles):
                            coefficients[i_bath,i_spin,i_pole] = -1j * temperature[i_bath][i_spin] * self.kappas_Fermi[i_pole]
                exponents=np.zeros((number_bath, number_spin, self.number_Pade_poles,2),complex)
                for i_bath in range(number_bath):
                    for i_spin in range(number_spin):
                        for i_pole in range(self.number_Pade_poles):
                            exponents[i_bath,i_spin,i_pole,0] = self.Xis_Fermi[i_pole]*temperature[i_bath][i_spin]-1j*chem_potentials[i_bath][i_spin]
                            exponents[i_bath,i_spin,i_pole,1] = self.Xis_Fermi[i_pole]*temperature[i_bath][i_spin]+1j*chem_potentials[i_bath][i_spin]
            else:
                coefficients=np.zeros((number_bath, number_spin, (self.number_Pade_poles+1)),complex)
                for i_bath in range(number_bath):
                    for i_spin in range(number_spin):
                        for i_pole in range(self.number_Pade_poles+1):
                            if i_pole==0:
                                #coefficients[i_bath,i_spin,i_pole] = np.pi *  broadenings[i_bath][i_spin] * self.Fermi_Pade_Approx(1j*broadenings[i_bath][i_spin],0.,temperature[i_bath][i_spin])
                                coefficients[i_bath,i_spin,i_pole] = broadenings[i_bath][i_spin] * self.Fermi_Pade_Approx(1j*broadenings[i_bath][i_spin],0.,temperature[i_bath][i_spin])/2. 
                            else:
                                #coefficients[i_bath,i_spin,i_pole] =  -1j * 2*np.pi * temperature[i_bath][i_spin] * self.kappas_Fermi[i_pole-1] * broadenings[i_bath][i_spin]**2 /( broadenings[i_bath][i_spin]**2 - (temperature[i_bath][i_spin]*self.Xis_Fermi[i_pole-1])**2 ) )
                                coefficients[i_bath,i_spin,i_pole] =  -1j * temperature[i_bath][i_spin] * self.kappas_Fermi[i_pole-1] * broadenings[i_bath][i_spin]**2 /( broadenings[i_bath][i_spin]**2 - (temperature[i_bath][i_spin]*self.Xis_Fermi[i_pole-1])**2 )  
                exponents=np.zeros((number_bath, number_spin, (self.number_Pade_poles+1),2),complex)
                for i_bath in range(number_bath):
                    for i_spin in range(number_spin):
                        for i_pole in range(self.number_Pade_poles+1):
                            if i_pole==0:
                                exponents[i_bath,i_spin,i_pole,0] = broadenings[i_bath][i_spin]-1j*chem_potentials[i_bath][i_spin]
                                exponents[i_bath,i_spin,i_pole,1] = broadenings[i_bath][i_spin]+1j*chem_potentials[i_bath][i_spin]
                            else:
                                exponents[i_bath,i_spin,i_pole,0] = self.Xis_Fermi[i_pole-1]*temperature[i_bath][i_spin]-1j*chem_potentials[i_bath][i_spin]
                                exponents[i_bath,i_spin,i_pole,1] = self.Xis_Fermi[i_pole-1]*temperature[i_bath][i_spin]+1j*chem_potentials[i_bath][i_spin]
            return coefficients, exponents 

if __name__=='__main__':
    number_pade_modes=5
    spec='Bose'
    pade=Pade_decomposition(number_pade_modes,spec)
    KAPPA,XI=pade.get_Pade_parameters()
    couplings=[1,1]
    broadenings=[10,10]
    Delta, coeffs, expos = pade.get_coefficients_exponents(1, couplings, broadenings, 2) 
    for i in range(len(KAPPA)):
        print(i,KAPPA[i],XI[i])
    if 'Fermi' in spec or 'fermi' in spec:
        print(pade.Fermi_V(energy=1.e-2,temperature=50.))
        print(pade.Fermi_Pade_Approx(energy=1.e-2,temperature=50.))
    elif 'Bose' in spec or 'bose' in spec:
        print(pade.Bose_V(energy=1.e-2,temperature=50.))
        print(pade.Bose_Pade_Approx(energy=1.e-2,temperature=50.))
        print(Delta)
        print(coeffs)
        print(expos)
        p_coeffs = np.sqrt( np.abs(coeffs) )
        m_coeffs = coeffs/p_coeffs
        print(p_coeffs)
        print(m_coeffs)
        p_coefficients = np.zeros( (2, number_pade_modes+1), dtype = complex)
        m_coefficients = np.zeros( (2, number_pade_modes+1), dtype = complex)
        for i in range(2):
            p_coefficients[i] = np.sqrt( np.abs(coeffs[i]) )
            m_coefficients[i] = coeffs[i]/p_coeffs[i]
        print(p_coefficients)
        print(m_coefficients)


    energies=np.arange(2.5,3.5,0.1)
    yy=[]
    for energy in energies:
        yy.append(pade.Fermi_V(energy,mu=3.0,temperature=300))

    import matplotlib.pyplot as plt
    plt.figure()
    plt.plot(energies,yy)
    plt.show()

