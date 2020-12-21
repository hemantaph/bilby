import numpy as np
import sympy as sp

C = 299792458.0

class Fn:
    def __init__(self, e, beta, iota):
        self.e = e
        self.beta = beta
        self.iota = iota
    

    def cplus(self, a, b):
        #(a,b)->(j,n)
        et = self.e
        beta = self.beta
        iota = self.iota
        Cp = np.zeros((10,3))

        Cp[0,-2] = et**3*( - (13/16)*sp.cos(2*beta) - 13/16*sp.cos(iota)**2*sp.cos(2*beta)) + et**5*( - (5/384)*sp.cos(2*beta) - 5/384*sp.cos(iota)**2*sp.cos(2*beta)) + et*(3/2*sp.cos(2*beta) + 3/2*sp.cos(iota)**2*sp.cos(2*beta))

        Cp[1,-2] = - 2*sp.cos(2*beta) - 2*sp.cos(iota)**2*sp.cos(2*beta) + et**4*( - (23/8)*sp.cos(2*beta) - 23/8*sp.cos(iota)**2*sp.cos(2*beta)) + et**6*(65/144*sp.cos(2*beta) + 65/144*sp.cos(iota)**2*sp.cos(2*beta)) + et**2*(5*sp.cos(2*beta) + 5*sp.cos(iota)**2*sp.cos(2*beta))

        Cp[2,-2] = et**5*( - (963/128)*sp.cos(2*beta) - 963/128*sp.cos(iota)**2*sp.cos(2*beta)) + et*( - (9/2)*sp.cos(2*beta) - 9/2*sp.cos(iota)**2*sp.cos(2*beta)) + et**3*(171/16*sp.cos(2*beta) + 171/16*sp.cos(iota)**2*sp.cos(2*beta))

        Cp[3,-2] = et**6*( - (101/6)*sp.cos(2*beta) - 101/6*sp.cos(iota)**2*sp.cos(2*beta)) + et**2*( - 8*sp.cos(2*beta) - 8*sp.cos(iota)**2*sp.cos(2*beta)) + et**4*(20*sp.cos(2*beta) + 20*sp.cos(iota)**2*sp.cos(2*beta))

        Cp[4,-2] = et**3*( - (625/48)*sp.cos(2*beta) - 625/48*sp.cos(iota)**2*sp.cos(2*beta)) + et**5*(26875/768*sp.cos(2*beta) + 26875/768*sp.cos(iota)**2*sp.cos(2*beta))

        Cp[5,-2] = et**4*( - (81/4)*sp.cos(2*beta) - 81/4*sp.cos(iota)**2*sp.cos(2*beta)) + et**6*(2349/40*sp.cos(2*beta) + 2349/40*sp.cos(iota)**2*sp.cos(2*beta))

        Cp[6,-2] = et**5*( - ((117649*sp.cos(2*beta))/3840) - (117649*sp.cos(iota)**2*sp.cos(2*beta))/3840)

        Cp[7,-2] = et**6*( - (2048/45)*sp.cos(2*beta) - 2048/45*sp.cos(iota)**2*sp.cos(2*beta))


        Cp[0,0] = et*sp.sin(iota)**2 - 1/8*et**3*sp.sin(iota)**2 + 1/192*et**5*sp.sin(iota)**2

        Cp[1,0] = et**2*sp.sin(iota)**2 - 1/3*et**4*sp.sin(iota)**2 + 1/24*et**6*sp.sin(iota)**2

        Cp[2,0] = 9/8*et**3*sp.sin(iota)**2 - 81/128*et**5*sp.sin(iota)**2

        Cp[3,0] = 4/3*et**4*sp.sin(iota)**2 - 16/15*et**6*sp.sin(iota)**2

        Cp[4,0] = 625/384*et**5*sp.sin(iota)**2

        Cp[5,0] = 81/40*et**6*sp.sin(iota)**2


        Cp[0,2] = et**5*(47/768*sp.cos(2*beta) + 47/768*sp.cos(iota)**2*sp.cos(2*beta)) + et**3*(7/48*sp.cos(2*beta) + 7/48*sp.cos(iota)**2*sp.cos(2*beta))

        Cp[1,2] = et**6*(11/240*sp.cos(2*beta) + 11/240*sp.cos(iota)**2*sp.cos(2*beta)) + et**4*(1/8*sp.cos(2*beta) + 1/8*sp.cos(iota)**2*sp.cos(2*beta))

        Cp[2,2] = et**5*((153*sp.cos(2*beta))/1280 + (153*sp.cos(iota)**2*sp.cos(2*beta))/1280)

        Cp[3,2] = et**6*(11/90*sp.cos(2*beta) + 11/90*sp.cos(iota)**2*sp.cos(2*beta))


        return(Cp[a-1,b])
    
    
    def splus(self, a, b):
        et = self.e
        beta = self.beta
        iota = self.iota
        Sp = np.zeros((10,3))

        Sp[0,-2] = et**3*( - (13/16)*sp.sin(2*beta) - 13/16*sp.cos(iota)**2*sp.sin(2*beta)) + et**5*( - (5/384)*sp.sin(2*beta) - 5/384*sp.cos(iota)**2*sp.sin(2*beta)) + et*(3/2*sp.sin(2*beta) + 3/2*sp.cos(iota)**2*sp.sin(2*beta))

        Sp[1,-2] = - 2*sp.sin(2*beta) - 2*sp.cos(iota)**2*sp.sin(2*beta) + et**4*( - (23/8)*sp.sin(2*beta) - 23/8*sp.cos(iota)**2*sp.sin(2*beta)) + et**6*(65/144*sp.sin(2*beta) + 65/144*sp.cos(iota)**2*sp.sin(2*beta)) + et**2*(5*sp.sin(2*beta) + 5*sp.cos(iota)**2*sp.sin(2*beta))

        Sp[2,-2] = et**5*( - (963/128)*sp.sin(2*beta) - 963/128*sp.cos(iota)**2*sp.sin(2*beta)) + et*( - (9/2)*sp.sin(2*beta) - 9/2*sp.cos(iota)**2*sp.sin(2*beta)) + et**3*(171/16*sp.sin(2*beta) + 171/16*sp.cos(iota)**2*sp.sin(2*beta))

        Sp[3,-2] = et**6*( - (101/6)*sp.sin(2*beta) - 101/6*sp.cos(iota)**2*sp.sin(2*beta)) + et**2*( - 8*sp.sin(2*beta) - 8*sp.cos(iota)**2*sp.sin(2*beta)) + et**4*(20*sp.sin(2*beta) + 20*sp.cos(iota)**2*sp.sin(2*beta))

        Sp[4,-2] = et**3*( - (625/48)*sp.sin(2*beta) - 625/48*sp.cos(iota)**2*sp.sin(2*beta)) + et**5*(26875/768*sp.sin(2*beta) + 26875/768*sp.cos(iota)**2*sp.sin(2*beta))

        Sp[5,-2] = et**4*( - (81/4)*sp.sin(2*beta) - 81/4*sp.cos(iota)**2*sp.sin(2*beta)) + et**6*(2349/40*sp.sin(2*beta) + 2349/40*sp.cos(iota)**2*sp.sin(2*beta))

        Sp[6,-2] = et**5*( - ((117649*sp.sin(2*beta))/3840) - (117649*sp.cos(iota)**2*sp.sin(2*beta))/3840)

        Sp[7,-2] = et**6*( - (2048/45)*sp.sin(2*beta) - 2048/45*sp.cos(iota)**2*sp.sin(2*beta))


        Sp[0,2] = et**3*( - (7/48)*sp.sin(2*beta) - 7/48*sp.cos(iota)**2*sp.sin(2*beta)) + et**5*( - (47/768)*sp.sin(2*beta) - 47/768*sp.cos(iota)**2*sp.sin(2*beta))

        Sp[1,2] = et**4*( - (1/8)*sp.sin(2*beta) - 1/8*sp.cos(iota)**2*sp.sin(2*beta)) + et**6*( - (11/240)*sp.sin(2*beta) - 11/240*sp.cos(iota)**2*sp.sin(2*beta))

        Sp[2,2] = et**5*( - ((153*sp.sin(2*beta))/1280) - (153*sp.cos(iota)**2*sp.sin(2*beta))/1280)

        Sp[3,2] = et**6*( - (11/90)*sp.sin(2*beta) - 11/90*sp.cos(iota)**2*sp.sin(2*beta))
        

        return(Sp[a-1,b])
    
    
    def ccross(self, a, b):
        et = self.e
        beta = self.beta
        iota = self.iota
        Cx = np.zeros((10,3))

        Cx[0,-2] = - 3*et*sp.cos(iota)*sp.sin(2*beta) + 13/8*et**3*sp.cos(iota)*sp.sin(2*beta) + 5/192*et**5*sp.cos(iota)*sp.sin(2*beta)

        Cx[1,-2] = 4*sp.cos(iota)*sp.sin(2*beta) - 10*et**2*sp.cos(iota)*sp.sin(2*beta) + 23/4*et**4*sp.cos(iota)*sp.sin(2*beta) - 65/72*et**6*sp.cos(iota)*sp.sin(2*beta)

        Cx[2,-2] = 9*et*sp.cos(iota)*sp.sin(2*beta) - 171/8*et**3*sp.cos(iota)*sp.sin(2*beta) + 963/64*et**5*sp.cos(iota)*sp.sin(2*beta)

        Cx[3,-2] = 16*et**2*sp.cos(iota)*sp.sin(2*beta) - 40*et**4*sp.cos(iota)*sp.sin(2*beta) + 101/3*et**6*sp.cos(iota)*sp.sin(2*beta)

        Cx[4,-2] = 625/24*et**3*sp.cos(iota)*sp.sin(2*beta) - 26875/384*et**5*sp.cos(iota)*sp.sin(2*beta)

        Cx[5,-2] = 81/2*et**4*sp.cos(iota)*sp.sin(2*beta) - 2349/20*et**6*sp.cos(iota)*sp.sin(2*beta)

        Cx[6,-2] = (117649*et**5*sp.cos(iota)*sp.sin(2*beta))/1920

        Cx[7,-2] = 4096/45*et**6*sp.cos(iota)*sp.sin(2*beta)


        Cx[0,2] = - (7/24)*et**3*sp.cos(iota)*sp.sin(2*beta) - 47/384*et**5*sp.cos(iota)*sp.sin(2*beta)

        Cx[1,2] = - (1/4)*et**4*sp.cos(iota)*sp.sin(2*beta) - 11/120*et**6*sp.cos(iota)*sp.sin(2*beta)

        Cx[2,2] = - (153/640)*et**5*sp.cos(iota)*sp.sin(2*beta)

        Cx[3,2] = - (11/45)*et**6*sp.cos(iota)*sp.sin(2*beta)
        

        return(Cx[a-1,b])
    
    
    def scross(self, a, b):
        et = self.e
        beta = self.beta
        iota = self.iota
        Sx = np.zeros((10,3))
        
        Sx[0,-2] = 3*et*sp.cos(iota)*sp.cos(2*beta) - 13/8*et**3*sp.cos(iota)*sp.cos(2*beta) - 5/192*et**5*sp.cos(iota)*sp.cos(2*beta)

        Sx[1,-2] = - 4*sp.cos(iota)*sp.cos(2*beta) + 10*et**2*sp.cos(iota)*sp.cos(2*beta) - 23/4*et**4*sp.cos(iota)*sp.cos(2*beta) + 65/72*et**6*sp.cos(iota)*sp.cos(2*beta)

        Sx[2,-2] = - 9*et*sp.cos(iota)*sp.cos(2*beta) + 171/8*et**3*sp.cos(iota)*sp.cos(2*beta) - 963/64*et**5*sp.cos(iota)*sp.cos(2*beta)

        Sx[3,-2] = - 16*et**2*sp.cos(iota)*sp.cos(2*beta) + 40*et**4*sp.cos(iota)*sp.cos(2*beta) - 101/3*et**6*sp.cos(iota)*sp.cos(2*beta)

        Sx[4,-2] = - (625/24)*et**3*sp.cos(iota)*sp.cos(2*beta) + 26875/384*et**5*sp.cos(iota)*sp.cos(2*beta)

        Sx[5,-2] = - (81/2)*et**4*sp.cos(iota)*sp.cos(2*beta) + 2349/20*et**6*sp.cos(iota)*sp.cos(2*beta)

        Sx[6,-2] = - ((117649*et**5*sp.cos(iota)*sp.cos(2*beta))/1920)

        Sx[7,-2] = - (4096/45)*et**6*sp.cos(iota)*sp.cos(2*beta)


        Sx[0,2] = - (7/24)*et**3*sp.cos(iota)*sp.cos(2*beta) - 47/384*et**5*sp.cos(iota)*sp.cos(2*beta)

        Sx[1,2] = - (1/4)*et**4*sp.cos(iota)*sp.cos(2*beta) - 11/120*et**6*sp.cos(iota)*sp.cos(2*beta)

        Sx[2,2] = - (153/640)*et**5*sp.cos(iota)*sp.cos(2*beta)

        Sx[3,2] = - (11/45)*et**6*sp.cos(iota)*sp.cos(2*beta)


        return(Sx[a-1,b])
    
    
    
    
    
    
    
    
    
    
    
    