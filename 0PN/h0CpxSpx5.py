import numpy as np
import sympy as sp

C = 299792458.0

class Fn:
    def __init__(self, e, beta, iota):
        self.e = e
        self.beta = beta
        self.iota = iota
    
    #it is plausible to take small sine reullt as 0
    #coz sine(pi) is 0, but np.sin(np.pi) is not 0
    #but u cannot the resolve the beta and iota angle between pi and np.pi 
    def sine_(self, x_):
        result_ = np.sin(float(x_)) 
        if abs(result_)<=4.898587196589413e-16:
            result_ = 0.0
            
        return(result_)
    
    def cosine_(self, x_):
        result_ = np.cos(float(x_)) 
        if abs(result_)<=1.8369701987210297e-16:
            result_ = 0.0
            
        return(result_)
    

    def cplus(self, a, b):
        #(a,b)->(j,n)
        sine = self.sine_
        cosine = self.cosine_
        
        et = self.e
        beta = self.beta
        iota = self.iota
        Cp = np.zeros((10,3))

        Cp[0,-2] = et**3*( - (13/16)*cosine(2*beta) - 13/16*cosine(iota)**2*cosine(2*beta)) + et**5*( - (5/384)*cosine(2*beta) - 5/384*cosine(iota)**2*cosine(2*beta)) + et*(3/2*cosine(2*beta) + 3/2*cosine(iota)**2*cosine(2*beta))

        Cp[1,-2] = - 2*cosine(2*beta) - 2*cosine(iota)**2*cosine(2*beta) + et**4*( - (23/8)*cosine(2*beta) - 23/8*cosine(iota)**2*cosine(2*beta)) + et**6*(65/144*cosine(2*beta) + 65/144*cosine(iota)**2*cosine(2*beta)) + et**2*(5*cosine(2*beta) + 5*cosine(iota)**2*cosine(2*beta))

        Cp[2,-2] = et**5*( - (963/128)*cosine(2*beta) - 963/128*cosine(iota)**2*cosine(2*beta)) + et*( - (9/2)*cosine(2*beta) - 9/2*cosine(iota)**2*cosine(2*beta)) + et**3*(171/16*cosine(2*beta) + 171/16*cosine(iota)**2*cosine(2*beta))

        Cp[3,-2] = et**6*( - (101/6)*cosine(2*beta) - 101/6*cosine(iota)**2*cosine(2*beta)) + et**2*( - 8*cosine(2*beta) - 8*cosine(iota)**2*cosine(2*beta)) + et**4*(20*cosine(2*beta) + 20*cosine(iota)**2*cosine(2*beta))

        Cp[4,-2] = et**3*( - (625/48)*cosine(2*beta) - 625/48*cosine(iota)**2*cosine(2*beta)) + et**5*(26875/768*cosine(2*beta) + 26875/768*cosine(iota)**2*cosine(2*beta))

        Cp[5,-2] = et**4*( - (81/4)*cosine(2*beta) - 81/4*cosine(iota)**2*cosine(2*beta)) + et**6*(2349/40*cosine(2*beta) + 2349/40*cosine(iota)**2*cosine(2*beta))

        Cp[6,-2] = et**5*( - ((117649*cosine(2*beta))/3840) - (117649*cosine(iota)**2*cosine(2*beta))/3840)

        Cp[7,-2] = et**6*( - (2048/45)*cosine(2*beta) - 2048/45*cosine(iota)**2*cosine(2*beta))


        Cp[0,0] = et*sine(iota)**2 - 1/8*et**3*sine(iota)**2 + 1/192*et**5*sine(iota)**2

        Cp[1,0] = et**2*sine(iota)**2 - 1/3*et**4*sine(iota)**2 + 1/24*et**6*sine(iota)**2

        Cp[2,0] = 9/8*et**3*sine(iota)**2 - 81/128*et**5*sine(iota)**2

        Cp[3,0] = 4/3*et**4*sine(iota)**2 - 16/15*et**6*sine(iota)**2

        Cp[4,0] = 625/384*et**5*sine(iota)**2

        Cp[5,0] = 81/40*et**6*sine(iota)**2


        Cp[0,2] = et**5*(47/768*cosine(2*beta) + 47/768*cosine(iota)**2*cosine(2*beta)) + et**3*(7/48*cosine(2*beta) + 7/48*cosine(iota)**2*cosine(2*beta))

        Cp[1,2] = et**6*(11/240*cosine(2*beta) + 11/240*cosine(iota)**2*cosine(2*beta)) + et**4*(1/8*cosine(2*beta) + 1/8*cosine(iota)**2*cosine(2*beta))

        Cp[2,2] = et**5*((153*cosine(2*beta))/1280 + (153*cosine(iota)**2*cosine(2*beta))/1280)

        Cp[3,2] = et**6*(11/90*cosine(2*beta) + 11/90*cosine(iota)**2*cosine(2*beta))


        return(Cp[a-1,b])
    
    
    def splus(self, a, b):
        sine = self.sine_
        cosine = self.cosine_
        
        et = self.e
        beta = self.beta
        iota = self.iota
        Sp = np.zeros((10,3))

        Sp[0,-2] = et**3*( - (13/16)*sine(2*beta) - 13/16*cosine(iota)**2*sine(2*beta)) + et**5*( - (5/384)*sine(2*beta) - 5/384*cosine(iota)**2*sine(2*beta)) + et*(3/2*sine(2*beta) + 3/2*cosine(iota)**2*sine(2*beta))

        Sp[1,-2] = - 2*sine(2*beta) - 2*cosine(iota)**2*sine(2*beta) + et**4*( - (23/8)*sine(2*beta) - 23/8*cosine(iota)**2*sine(2*beta)) + et**6*(65/144*sine(2*beta) + 65/144*cosine(iota)**2*sine(2*beta)) + et**2*(5*sine(2*beta) + 5*cosine(iota)**2*sine(2*beta))

        Sp[2,-2] = et**5*( - (963/128)*sine(2*beta) - 963/128*cosine(iota)**2*sine(2*beta)) + et*( - (9/2)*sine(2*beta) - 9/2*cosine(iota)**2*sine(2*beta)) + et**3*(171/16*sine(2*beta) + 171/16*cosine(iota)**2*sine(2*beta))

        Sp[3,-2] = et**6*( - (101/6)*sine(2*beta) - 101/6*cosine(iota)**2*sine(2*beta)) + et**2*( - 8*sine(2*beta) - 8*cosine(iota)**2*sine(2*beta)) + et**4*(20*sine(2*beta) + 20*cosine(iota)**2*sine(2*beta))

        Sp[4,-2] = et**3*( - (625/48)*sine(2*beta) - 625/48*cosine(iota)**2*sine(2*beta)) + et**5*(26875/768*sine(2*beta) + 26875/768*cosine(iota)**2*sine(2*beta))

        Sp[5,-2] = et**4*( - (81/4)*sine(2*beta) - 81/4*cosine(iota)**2*sine(2*beta)) + et**6*(2349/40*sine(2*beta) + 2349/40*cosine(iota)**2*sine(2*beta))

        Sp[6,-2] = et**5*( - ((117649*sine(2*beta))/3840) - (117649*cosine(iota)**2*sine(2*beta))/3840)

        Sp[7,-2] = et**6*( - (2048/45)*sine(2*beta) - 2048/45*cosine(iota)**2*sine(2*beta))


        Sp[0,2] = et**3*( - (7/48)*sine(2*beta) - 7/48*cosine(iota)**2*sine(2*beta)) + et**5*( - (47/768)*sine(2*beta) - 47/768*cosine(iota)**2*sine(2*beta))

        Sp[1,2] = et**4*( - (1/8)*sine(2*beta) - 1/8*cosine(iota)**2*sine(2*beta)) + et**6*( - (11/240)*sine(2*beta) - 11/240*cosine(iota)**2*sine(2*beta))

        Sp[2,2] = et**5*( - ((153*sine(2*beta))/1280) - (153*cosine(iota)**2*sine(2*beta))/1280)

        Sp[3,2] = et**6*( - (11/90)*sine(2*beta) - 11/90*cosine(iota)**2*sine(2*beta))
        

        return(Sp[a-1,b])
    
    
    def ccross(self, a, b):
        sine = self.sine_
        cosine = self.cosine_
        
        et = self.e
        beta = self.beta
        iota = self.iota
        Cx = np.zeros((10,3))

        Cx[0,-2] = - 3*et*cosine(iota)*sine(2*beta) + 13/8*et**3*cosine(iota)*sine(2*beta) + 5/192*et**5*cosine(iota)*sine(2*beta)

        Cx[1,-2] = 4*cosine(iota)*sine(2*beta) - 10*et**2*cosine(iota)*sine(2*beta) + 23/4*et**4*cosine(iota)*sine(2*beta) - 65/72*et**6*cosine(iota)*sine(2*beta)

        Cx[2,-2] = 9*et*cosine(iota)*sine(2*beta) - 171/8*et**3*cosine(iota)*sine(2*beta) + 963/64*et**5*cosine(iota)*sine(2*beta)

        Cx[3,-2] = 16*et**2*cosine(iota)*sine(2*beta) - 40*et**4*cosine(iota)*sine(2*beta) + 101/3*et**6*cosine(iota)*sine(2*beta)

        Cx[4,-2] = 625/24*et**3*cosine(iota)*sine(2*beta) - 26875/384*et**5*cosine(iota)*sine(2*beta)

        Cx[5,-2] = 81/2*et**4*cosine(iota)*sine(2*beta) - 2349/20*et**6*cosine(iota)*sine(2*beta)

        Cx[6,-2] = (117649*et**5*cosine(iota)*sine(2*beta))/1920

        Cx[7,-2] = 4096/45*et**6*cosine(iota)*sine(2*beta)


        Cx[0,2] = - (7/24)*et**3*cosine(iota)*sine(2*beta) - 47/384*et**5*cosine(iota)*sine(2*beta)

        Cx[1,2] = - (1/4)*et**4*cosine(iota)*sine(2*beta) - 11/120*et**6*cosine(iota)*sine(2*beta)

        Cx[2,2] = - (153/640)*et**5*cosine(iota)*sine(2*beta)

        Cx[3,2] = - (11/45)*et**6*cosine(iota)*sine(2*beta)
        

        return(Cx[a-1,b])
    
    
    def scross(self, a, b):
        sine = self.sine_
        cosine = self.cosine_
        
        et = self.e
        beta = self.beta
        iota = self.iota
        Sx = np.zeros((10,3))
        
        Sx[0,-2] = 3*et*cosine(iota)*cosine(2*beta) - 13/8*et**3*cosine(iota)*cosine(2*beta) - 5/192*et**5*cosine(iota)*cosine(2*beta)

        Sx[1,-2] = - 4*cosine(iota)*cosine(2*beta) + 10*et**2*cosine(iota)*cosine(2*beta) - 23/4*et**4*cosine(iota)*cosine(2*beta) + 65/72*et**6*cosine(iota)*cosine(2*beta)

        Sx[2,-2] = - 9*et*cosine(iota)*cosine(2*beta) + 171/8*et**3*cosine(iota)*cosine(2*beta) - 963/64*et**5*cosine(iota)*cosine(2*beta)

        Sx[3,-2] = - 16*et**2*cosine(iota)*cosine(2*beta) + 40*et**4*cosine(iota)*cosine(2*beta) - 101/3*et**6*cosine(iota)*cosine(2*beta)

        Sx[4,-2] = - (625/24)*et**3*cosine(iota)*cosine(2*beta) + 26875/384*et**5*cosine(iota)*cosine(2*beta)

        Sx[5,-2] = - (81/2)*et**4*cosine(iota)*cosine(2*beta) + 2349/20*et**6*cosine(iota)*cosine(2*beta)

        Sx[6,-2] = - ((117649*et**5*cosine(iota)*cosine(2*beta))/1920)

        Sx[7,-2] = - (4096/45)*et**6*cosine(iota)*cosine(2*beta)


        Sx[0,2] = - (7/24)*et**3*cosine(iota)*cosine(2*beta) - 47/384*et**5*cosine(iota)*cosine(2*beta)

        Sx[1,2] = - (1/4)*et**4*cosine(iota)*cosine(2*beta) - 11/120*et**6*cosine(iota)*cosine(2*beta)

        Sx[2,2] = - (153/640)*et**5*cosine(iota)*cosine(2*beta)

        Sx[3,2] = - (11/45)*et**6*cosine(iota)*cosine(2*beta)


        return(Sx[a-1,b])
    
    
    
    
    
    
    
    
    
    
    
    