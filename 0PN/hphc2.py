import numpy as np
import matplotlib.pyplot as plt

import h0CpxSpx

C = 299792458.
G = 6.67408*1e-11
Mo = 1.989*1e30
Mpc = 3.086*1e22

class Fn:
    def __init__(self, iota_, beta_, D_, m1_, m2_, f_, f0_, Fp_, Fc_, et0_, phic_, tc_ ):

        self.iota_ = iota_
        self.beta_ = beta_
        self.D_ = D_*Mpc
        self.m1_ = m1_*Mo
        self.m2_ = m2_*Mo
        self.f0_ = f0_
        self.Fp_ = Fp_
        self.Fc_ = Fc_
        self.et0_ = et0_
        self.phic_ = phic_
        self.tc_ = tc_
        self.f_ = f_

    #defining unit-step function 
    def unitstep(self,lp,ff,f):
        if lp*ff-2*f>=0:
            return(1)
        else:
            return(0)
    
    #check which of the harmonics at a particular frequency (2*Pi*Omega=f/l) exceeds the  
    def uniarray(self, l, ff, f):
        lx = np.arange(1,11,1)
        for i in range(10):
            lx[i] = self.unitstep(l[i], ff, f)

        return(lx)
    
    #eccentricity
    def eccn(self):
        f = self.f_
        f0 = self.f0_
        et0 = self.et0_
        chi = f/f0
        
        et = et0**3*( - (3323/(1824*chi**(19/6))) + 3323/(1824*chi**(19/18)) ) + \
        et0**5*( 50259743/(6653952*chi**(95/18)) - 11042329/(1108992*chi**(19/6)) + 15994231/(6653952*chi**(19/18)) ) + \
        et0/chi**(19/18)
            
        return(et)
    
        #fourier phase
    def psi_l(self):
        l = np.arange(1,11,1)
        phic = self.phic_
        tc = self.tc_ 
        m1 = self.m1_
        m2 = self.m2_
        M = m1 + m2
        eta = (m1*m2)/(M**2)
        f = self.f_
        f0 = self.f0_
        et0 = self.et0_
        chi = f/f0       

        psi = l*phic - 2*np.pi*f*tc - 3/(128*eta) * ((G*M*np.pi*f)/C**3)**(-5/3) * \
        (l/2)**(8/3) *( 1 + et0**6 *( -(75356125/(3326976 *chi**(19/3)))+\
                                     (17355248095/(455518464 *chi**(38/9)))-\
                                     (1326481225/(101334144 *chi**(19/9))) ) + \
                       et0**4 *( (5222765/(998944 *chi**(38/9))) - (2608555/(444448 *chi**(19/9))) ) - \
                       (2355 *et0**2)/(1462 *chi**(19/9)) ) 
        
        return(psi)
    
    def xi_l(self):
        et = self.eccn()
        iota = self.iota_
        beta = self.beta_
        Fp = self.Fp_
        Fc = self.Fc_
        
         #calling class for Cx C+ Sx S+
        fn = h0CpxSpx.Fn(et,iota,beta)

        #to find Xi
        Gamma_l = Fp*fn.cplus() + Fc*fn.ccross()
        Sigma_l = Fp*fn.splus() + Fc*fn.scross()

        al = np.sign(Gamma_l)*np.sqrt(Gamma_l**2 + Sigma_l**2)

        phil = np.zeros((10,3))
        for i in range(10):
            for j in range(3):
                if Gamma_l[i][j]==0:
                    phil[i][j] = -np.sign(Sigma_l[i][j])*np.pi/2 
                else:
                    phil[i][j] = np.arctan(- (Sigma_l[i][j]/Gamma_l[i][j]))            
        
        numerator = (1-et**2)**(7/4)
        denomitor = ( 1 + (73/24)*et**2 + (37/96)*et**4 )**(1/2)
        xil = (numerator/denomitor)*al*np.exp(-1j*phil)    
        
        return(xil)
 
    def htilde0(self):
        
        k = 
        
        xil = self.xi_l0()
        psi = self.psi_l()
        l = np.arange(1,11,1)
        f = self.f_
        m1 = self.m1_
        m2 = self.m2_
        M = m1 + m2
        D = self.D_
        ff = C**3/(G*M*np.pi*6**(3/2))
        eta = (m1*m2)/(M**2)
        
        #frequency domain waveform
        s = np.sum( xil[:,0]*((l/2)**(2/3))*np.exp( -1j*(np.pi/4 + psi(l,0)) )*self.uniarray( (l-(l+0)*(k/(1+k))) ,ff,f) ) + \
        np.sum( xil[:,-2]*((l/2)**(2/3))*np.exp( -1j*(np.pi/4 + psi(l,-2)) )*self.uniarray( (l-(l-2)*(k/(1+k))) ,ff,f) ) + \
        np.sum( xil[:,2]*((l/2)**(2/3))*np.exp( -1j*(np.pi/4 + psi(l,2)) )*self.uniarray( (l-(l+2)*(k/(1+k))) ,ff,f) )
        
        hf0 = ((5*np.pi*eta)/384)**(1/2) * (G**2*M**2)/(C**5*D)*( ((G*M*np.pi*f)/C**3)**(-7/6) )*np.sum(s)
        
        return(hf0)
    
    
'''   
    def h_tilde(self):

        iota = self.iota_
        beta = self.beta_
        D = self.D_
        m1 = self.m1_
        m2 = self.m2_
        f = self.f_
        f0 = self.f0_
        Fp = self.Fp_
        Fc = self.Fc_
        et0 = self.et0_
        phic = self.phic_
        tc = self.tc_        
        
        M = m1 + m2
        eta = (m1*m2)/(M**2)
        chi = f/f0
        ff = C**3/(G*M*np.pi*6**(3/2))
        l = np.arange(1,11,1)

        #eccentricity
        et = et0**3*( - (3323/(1824*chi**(19/6))) + 3323/(1824*chi**(19/18)) ) + \
        et0**5*( 50259743/(6653952*chi**(95/18)) - 11042329/(1108992*chi**(19/6)) + 15994231/(6653952*chi**(19/18)) ) + et0/chi**(19/18)

        #calling class for Cx C+ Sx S+
        fn = h0CpxSpx.Fn(et,iota,beta)

        #to find Xi
        Gamma_l = Fp*fn.cplus() + Fc*fn.ccross()
        Sigma_l = Fp*fn.splus() + Fc*fn.scross()

        al = np.sign(Gamma_l)*np.sqrt(Gamma_l**2 + Sigma_l**2)

        phil = np.zeros((10,3))
        for i in range(10):
            for j in range(3):
                if Gamma_l[i][j]==0:
                    phil[i][j] = -np.sign(Sigma_l[i][j])*np.pi/2 
                else:
                    phil[i][j] = np.arctan(- (Sigma_l[i][j]/Gamma_l[i][j]))

        numerator = (1-et**2)**(7/4)
        denomitor = ( 1 + (73/24)*et**2 + (37/96)*et**4 )**(1/2)
        xil = (numerator/denomitor)*al*np.exp(-1j*phil)                

        #fourier phase
        psi = l*phic - 2*np.pi*f*tc - 3/(128*eta) * ((G*M*np.pi*f)/C**3)**(-5/3) * \
        (l/2)**(8/3) *( 1 + et0**6 *( -(75356125/(3326976 *chi**(19/3)))+\
                                     (17355248095/(455518464 *chi**(38/9)))-\
                                     (1326481225/(101334144 *chi**(19/9))) ) + \
                       et0**4 *( (5222765/(998944 *chi**(38/9))) - (2608555/(444448 *chi**(19/9))) ) - \
                       (2355 *et0**2)/(1462 *chi**(19/9)) ) 


        #frequency domain waveform
        s = np.sum( xil[:,0]*((l/2)**(2/3))*np.exp( -1j*(np.pi/4 + psi) )*self.uniarray(l,ff,f) ) + \
        np.sum( xil[:,-2]*((l/2)**(2/3))*np.exp( -1j*(np.pi/4 + psi) )*self.uniarray(l,ff,f) ) + \
        np.sum( xil[:,2]*((l/2)**(2/3))*np.exp( -1j*(np.pi/4 + psi) )*self.uniarray(l,ff,f) )
        hf = ((5*np.pi*eta)/384)**(1/2) * (G**2*M**2)/(C**5*D)*( ((G*M*np.pi*f)/C**3)**(-7/6) )*np.sum(s)
        
        return(hf)
'''        
        
        