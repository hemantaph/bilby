import numpy as np
import sympy as sp

C = 299792458.0

class Fn:
    def __init__(self, e, beta, iota):      
        self.e = e
        self.beta = beta
        self.iota = iota
    
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
        sine = self.sine_
        cosine = self.cosine_
        
        e = self.e
        beta = self.beta
        iota = self.iota
        Cp = np.zeros((10,5))

        Cp[0,1]=e**6*((109*cosine(beta)*sine(iota))/12288+(265*cosine(iota)**2*cosine(beta)*sine(iota))/12288)+e**4*(1/64*cosine(beta)*sine(iota)+7/192*cosine(iota)**2*cosine(beta)*sine(iota))+e**2*(-(9/32)*cosine(beta)*sine(iota)+11/32*cosine(iota)**2*cosine(beta)*sine(iota))
        Cp[1,1]=e**5*(1/16*cosine(beta)*sine(iota)-1/48*cosine(iota)**2*cosine(beta)*sine(iota))+e**3*(-(1/3)*cosine(beta)*sine(iota)+1/3*cosine(iota)**2*cosine(beta)*sine(iota))
        Cp[2,1]=e**6*((369*cosine(beta)*sine(iota))/2560-(243*cosine(iota)**2*cosine(beta)*sine(iota))/2560)+e**4*(-(207/512)*cosine(beta)*sine(iota)+189/512*cosine(iota)**2*cosine(beta)*sine(iota))
        Cp[3,1]=e**5*(-(1/2)*cosine(beta)*sine(iota)+13/30*cosine(iota)**2*cosine(beta)*sine(iota))
        Cp[4,1]=e**6*(-((23125*cosine(beta)*sine(iota))/36864)+(19375*cosine(iota)**2*cosine(beta)*sine(iota))/36864)


        Cp[0,-1]=-(5/4)*cosine(beta)*sine(iota)-1/4*cosine(iota)**2*cosine(beta)*sine(iota)+e**2*(3/2*cosine(beta)*sine(iota)-1/2*cosine(iota)**2*cosine(beta)*sine(iota))+e**6*((41*cosine(beta)*sine(iota))/2304+(37*cosine(iota)**2*cosine(beta)*sine(iota))/2304)+e**4*(-(51/256)*cosine(beta)*sine(iota)+41/256*cosine(iota)**2*cosine(beta)*sine(iota))
        Cp[1,-1]=e**3*(4*cosine(beta)*sine(iota)-2*cosine(iota)**2*cosine(beta)*sine(iota))+e**5*(-(19/16)*cosine(beta)*sine(iota)+35/48*cosine(iota)**2*cosine(beta)*sine(iota))+e*(-3*cosine(beta)*sine(iota)+cosine(iota)**2*cosine(beta)*sine(iota))
        Cp[2,-1]=e**4*(531/64*cosine(beta)*sine(iota)-297/64*cosine(iota)**2*cosine(beta)*sine(iota))+e**6*(-((15399*cosine(beta)*sine(iota))/4096)+(9477*cosine(iota)**2*cosine(beta)*sine(iota))/4096)+e**2*(-(171/32)*cosine(beta)*sine(iota)+81/32*cosine(iota)**2*cosine(beta)*sine(iota))
        Cp[3,-1]=e**5*(31/2*cosine(beta)*sine(iota)-55/6*cosine(iota)**2*cosine(beta)*sine(iota))+e**3*(-(26/3)*cosine(beta)*sine(iota)+14/3*cosine(iota)**2*cosine(beta)*sine(iota))
        Cp[4,-1]=e**6*((41875*cosine(beta)*sine(iota))/1536-(25625*cosine(iota)**2*cosine(beta)*sine(iota))/1536)+e**4*(-(6875/512)*cosine(beta)*sine(iota)+(11875*cosine(iota)**2*cosine(beta)*sine(iota))/1536)
        Cp[5,-1]=e**5*(-(81/4)*cosine(beta)*sine(iota)+243/20*cosine(iota)**2*cosine(beta)*sine(iota))
        Cp[6,-1]=e**6*(-((5529503*cosine(beta)*sine(iota))/184320)+(3411821*cosine(iota)**2*cosine(beta)*sine(iota))/184320)
        
        
        Cp[0,3]=e**4*(-(25/512)*cosine(3*beta)*sine(iota)-25/512*cosine(iota)**2*cosine(3*beta)*sine(iota))+e**6*(-((179*cosine(3*beta)*sine(iota))/7680)-(179*cosine(iota)**2*cosine(3*beta)*sine(iota))/7680)
        Cp[1,3]=e**5*(-(13/240)*cosine(3*beta)*sine(iota)-13/240*cosine(iota)**2*cosine(3*beta)*sine(iota))
        Cp[2,3]=e**6*(-((1233*cosine(3*beta)*sine(iota))/20480)-(1233*cosine(iota)**2*cosine(3*beta)*sine(iota))/20480)
      
    
        Cp[0,-3]=e**4*(-(65/192)*cosine(3*beta)*sine(iota)-65/192*cosine(iota)**2*cosine(3*beta)*sine(iota))+e**6*(-((19*cosine(3*beta)*sine(iota))/12288)-(19*cosine(iota)**2*cosine(3*beta)*sine(iota))/12288)+e**2*(19/32*cosine(3*beta)*sine(iota)+19/32*cosine(iota)**2*cosine(3*beta)*sine(iota))
        Cp[1,-3]=e*(-3*cosine(3*beta)*sine(iota)-3*cosine(iota)**2*cosine(3*beta)*sine(iota))+e**5*(-(133/48)*cosine(3*beta)*sine(iota)-133/48*cosine(iota)**2*cosine(3*beta)*sine(iota))+e**3*(11/2*cosine(3*beta)*sine(iota)+11/2*cosine(iota)**2*cosine(3*beta)*sine(iota))
        Cp[2,-3]=9/4*cosine(3*beta)*sine(iota)+9/4*cosine(iota)**2*cosine(3*beta)*sine(iota)+e**2*(-(27/2)*cosine(3*beta)*sine(iota)-27/2*cosine(iota)**2*cosine(3*beta)*sine(iota))+e**6*(-(3141/256)*cosine(3*beta)*sine(iota)-3141/256*cosine(iota)**2*cosine(3*beta)*sine(iota))+e**4*(5319/256*cosine(3*beta)*sine(iota)+5319/256*cosine(iota)**2*cosine(3*beta)*sine(iota))
        Cp[3,-3]=e**3*(-38*cosine(3*beta)*sine(iota)-38*cosine(iota)**2*cosine(3*beta)*sine(iota))+e*(8*cosine(3*beta)*sine(iota)+8*cosine(iota)**2*cosine(3*beta)*sine(iota))+e**5*(176/3*cosine(3*beta)*sine(iota)+176/3*cosine(iota)**2*cosine(3*beta)*sine(iota))
        Cp[4,-3]=e**4*(-(5625/64)*cosine(3*beta)*sine(iota)-5625/64*cosine(iota)**2*cosine(3*beta)*sine(iota))+e**2*(625/32*cosine(3*beta)*sine(iota)+625/32*cosine(iota)**2*cosine(3*beta)*sine(iota))+e**6*((1746875*cosine(3*beta)*sine(iota))/12288+(1746875*cosine(iota)**2*cosine(3*beta)*sine(iota))/12288)
        Cp[5,-3]=e**5*(-(729/4)*cosine(3*beta)*sine(iota)-729/4*cosine(iota)**2*cosine(3*beta)*sine(iota))+e**3*(81/2*cosine(3*beta)*sine(iota)+81/2*cosine(iota)**2*cosine(3*beta)*sine(iota))
        Cp[6,-3]=e**6*(-((2705927*cosine(3*beta)*sine(iota))/7680)-(2705927*cosine(iota)**2*cosine(3*beta)*sine(iota))/7680)+e**4*((117649*cosine(3*beta)*sine(iota))/1536+(117649*cosine(iota)**2*cosine(3*beta)*sine(iota))/1536)
        Cp[7,-3]=e**5*(2048/15*cosine(3*beta)*sine(iota)+2048/15*cosine(iota)**2*cosine(3*beta)*sine(iota))
        Cp[8,-3]=e**6*((4782969*cosine(3*beta)*sine(iota))/20480+(4782969*cosine(iota)**2*cosine(3*beta)*sine(iota))/20480)  


        return(Cp[a-1,b])
    
    
    def splus(self, a, b):
        sine = self.sine_
        cosine = self.cosine_
        
        e = self.e
        beta = self.beta
        iota = self.iota
        Sp = np.zeros((10,5))

        Sp[0,-1] =  - (5/4)*sine(iota)*sine(beta) - 1/4*cosine(iota)**2*sine(iota)*sine(beta) + e**2*(3/2*sine(iota)*sine(beta) - 1/2*cosine(iota)**2*sine(iota)*sine(beta)) + e**6*((41*sine(iota)*sine(beta))/2304 + (37*cosine(iota)**2*sine(iota)*sine(beta))/2304) + e**4*( - (51/256)*sine(iota)*sine(beta) + 41/256*cosine(iota)**2*sine(iota)*sine(beta))

        Sp[1,-1] = e**3*(4*sine(iota)*sine(beta) - 2*cosine(iota)**2*sine(iota)*sine(beta)) + e**5*( - (19/16)*sine(iota)*sine(beta) + 35/48*cosine(iota)**2*sine(iota)*sine(beta)) + e*( - 3*sine(iota)*sine(beta) + cosine(iota)**2*sine(iota)*sine(beta))

        Sp[2,-1] = e**4*(531/64*sine(iota)*sine(beta) - 297/64*cosine(iota)**2*sine(iota)*sine(beta)) + e**6*( - ((15399*sine(iota)*sine(beta))/4096) + (9477*cosine(iota)**2*sine(iota)*sine(beta))/4096) + e**2*( - (171/32)*sine(iota)*sine(beta) + 81/32*cosine(iota)**2*sine(iota)*sine(beta))

        Sp[3,-1] = e**5*(31/2*sine(iota)*sine(beta) - 55/6*cosine(iota)**2*sine(iota)*sine(beta)) + e**3*( - (26/3)*sine(iota)*sine(beta) + 14/3*cosine(iota)**2*sine(iota)*sine(beta))

        Sp[4,-1] = e**6*((41875*sine(iota)*sine(beta))/1536 - (25625*cosine(iota)**2*sine(iota)*sine(beta))/1536) + e**4*( - (6875/512)*sine(iota)*sine(beta) + (11875*cosine(iota)**2*sine(iota)*sine(beta))/1536)

        Sp[5,-1] = e**5*( - (81/4)*sine(iota)*sine(beta) + 243/20*cosine(iota)**2*sine(iota)*sine(beta))

        Sp[6,-1] = e**6*( - ((5529503*sine(iota)*sine(beta))/184320) + (3411821*cosine(iota)**2*sine(iota)*sine(beta))/184320)


        Sp[0,1] = e**2*(9/32*sine(iota)*sine(beta) - 11/32*cosine(iota)**2*sine(iota)*sine(beta)) + e**4*( - (1/64)*sine(iota)*sine(beta) - 7/192*cosine(iota)**2*sine(iota)*sine(beta)) + e**6*( - ((109*sine(iota)*sine(beta))/12288) - (265*cosine(iota)**2*sine(iota)*sine(beta))/12288)

        Sp[1,1] = e**3*(1/3*sine(iota)*sine(beta) - 1/3*cosine(iota)**2*sine(iota)*sine(beta)) + e**5*( - (1/16)*sine(iota)*sine(beta) + 1/48*cosine(iota)**2*sine(iota)*sine(beta))

        Sp[2,1] = e**4*(207/512*sine(iota)*sine(beta) - 189/512*cosine(iota)**2*sine(iota)*sine(beta)) + e**6*( - ((369*sine(iota)*sine(beta))/2560) + (243*cosine(iota)**2*sine(iota)*sine(beta))/2560)

        Sp[3,1] = e**5*(1/2*sine(iota)*sine(beta) - 13/30*cosine(iota)**2*sine(iota)*sine(beta))

        Sp[4,1] = e**6*((23125*sine(iota)*sine(beta))/36864 - (19375*cosine(iota)**2*sine(iota)*sine(beta))/36864)


        Sp[0,-3] = e**4*( - (65/192)*sine(iota)*sine(3*beta) - 65/192*cosine(iota)**2*sine(iota)*sine(3*beta)) + e**6*( - ((19*sine(iota)*sine(3*beta))/12288) - (19*cosine(iota)**2*sine(iota)*sine(3*beta))/12288) + e**2*(19/32*sine(iota)*sine(3*beta) + 19/32*cosine(iota)**2*sine(iota)*sine(3*beta))

        Sp[1,-3] = e*( - 3*sine(iota)*sine(3*beta) - 3*cosine(iota)**2*sine(iota)*sine(3*beta)) + e**5*( - (133/48)*sine(iota)*sine(3*beta) - 133/48*cosine(iota)**2*sine(iota)*sine(3*beta)) + e**3*(11/2*sine(iota)*sine(3*beta) + 11/2*cosine(iota)**2*sine(iota)*sine(3*beta))

        Sp[2,-3] = 9/4*sine(iota)*sine(3*beta) + 9/4*cosine(iota)**2*sine(iota)*sine(3*beta) + e**2*( - (27/2)*sine(iota)*sine(3*beta) - 27/2*cosine(iota)**2*sine(iota)*sine(3*beta)) + e**6*( - (3141/256)*sine(iota)*sine(3*beta) - 3141/256*cosine(iota)**2*sine(iota)*sine(3*beta)) + e**4*(5319/256*sine(iota)*sine(3*beta) + 5319/256*cosine(iota)**2*sine(iota)*sine(3*beta))

        Sp[3,-3] = e**3*( - 38*sine(iota)*sine(3*beta) - 38*cosine(iota)**2*sine(iota)*sine(3*beta)) + e*(8*sine(iota)*sine(3*beta) + 8*cosine(iota)**2*sine(iota)*sine(3*beta)) + e**5*(176/3*sine(iota)*sine(3*beta) + 176/3*cosine(iota)**2*sine(iota)*sine(3*beta))

        Sp[4,-3] = e**4*( - (5625/64)*sine(iota)*sine(3*beta) - 5625/64*cosine(iota)**2*sine(iota)*sine(3*beta)) + e**2*(625/32*sine(iota)*sine(3*beta) + 625/32*cosine(iota)**2*sine(iota)*sine(3*beta)) + e**6*((1746875*sine(iota)*sine(3*beta))/12288 + (1746875*cosine(iota)**2*sine(iota)*sine(3*beta))/12288)

        Sp[5,-3] = e**5*( - (729/4)*sine(iota)*sine(3*beta) - 729/4*cosine(iota)**2*sine(iota)*sine(3*beta)) + e**3*(81/2*sine(iota)*sine(3*beta) + 81/2*cosine(iota)**2*sine(iota)*sine(3*beta))

        Sp[6,-3] = e**6*( - ((2705927*sine(iota)*sine(3*beta))/7680) - (2705927*cosine(iota)**2*sine(iota)*sine(3*beta))/7680) + e**4*((117649*sine(iota)*sine(3*beta))/1536 + (117649*cosine(iota)**2*sine(iota)*sine(3*beta))/1536)

        Sp[7,-3] = e**5*(2048/15*sine(iota)*sine(3*beta) + 2048/15*cosine(iota)**2*sine(iota)*sine(3*beta))

        Sp[8,-3] = e**6*((4782969*sine(iota)*sine(3*beta))/20480 + (4782969*cosine(iota)**2*sine(iota)*sine(3*beta))/20480)


        Sp[0,3] = e**6*((179*sine(iota)*sine(3*beta))/7680 + (179*cosine(iota)**2*sine(iota)*sine(3*beta))/7680) + e**4*(25/512*sine(iota)*sine(3*beta) + 25/512*cosine(iota)**2*sine(iota)*sine(3*beta))

        Sp[1,3] = e**5*(13/240*sine(iota)*sine(3*beta) + 13/240*cosine(iota)**2*sine(iota)*sine(3*beta))

        Sp[2,3] = e**6*((1233*sine(iota)*sine(3*beta))/20480 + (1233*cosine(iota)**2*sine(iota)*sine(3*beta))/20480)

        return(Sp[a-1,b])
    
    
    def ccross(self, a, b):
        sine = self.sine_
        cosine = self.cosine_
        
        e = self.e
        beta = self.beta
        iota = self.iota
        Cx = np.zeros((10,5))

        Cx[0,1]=-(1/16)*e**2*cosine(iota)*sine(iota)*sine(beta)-5/96*e**4*cosine(iota)*sine(iota)*sine(beta)-(187*e**6*cosine(iota)*sine(iota)*sine(beta))/6144
        Cx[1,1]=-(1/24)*e**5*cosine(iota)*sine(iota)*sine(beta)
        Cx[2,1]=9/256*e**4*cosine(iota)*sine(iota)*sine(beta)-(63*e**6*cosine(iota)*sine(iota)*sine(beta))/1280
        Cx[3,1]=1/15*e**5*cosine(iota)*sine(iota)*sine(beta)
        Cx[4,1]=(625*e**6*cosine(iota)*sine(iota)*sine(beta))/6144
        
        
        Cx[0,-1]=3/2*cosine(iota)*sine(iota)*sine(beta)-e**2*cosine(iota)*sine(iota)*sine(beta)+5/128*e**4*cosine(iota)*sine(iota)*sine(beta)-13/384*e**6*cosine(iota)*sine(iota)*sine(beta)
        Cx[1,-1]=2*e*cosine(iota)*sine(iota)*sine(beta)-2*e**3*cosine(iota)*sine(iota)*sine(beta)+11/24*e**5*cosine(iota)*sine(iota)*sine(beta)
        Cx[2,-1]=45/16*e**2*cosine(iota)*sine(iota)*sine(beta)-117/32*e**4*cosine(iota)*sine(iota)*sine(beta)+(2961*e**6*cosine(iota)*sine(iota)*sine(beta))/2048
        Cx[3,-1]=4*e**3*cosine(iota)*sine(iota)*sine(beta)-19/3*e**5*cosine(iota)*sine(iota)*sine(beta)
        Cx[4,-1]=4375/768*e**4*cosine(iota)*sine(iota)*sine(beta)-8125/768*e**6*cosine(iota)*sine(iota)*sine(beta)
        Cx[5,-1]=81/10*e**5*cosine(iota)*sine(iota)*sine(beta)
        Cx[6,-1]=(117649*e**6*cosine(iota)*sine(iota)*sine(beta))/10240 
        
        
        Cx[0,3]=25/256*e**4*cosine(iota)*sine(iota)*sine(3*beta)+(179*e**6*cosine(iota)*sine(iota)*sine(3*beta))/3840
        Cx[1,3]=13/120*e**5*cosine(iota)*sine(iota)*sine(3*beta)
        Cx[2,3]=(1233*e**6*cosine(iota)*sine(iota)*sine(3*beta))/10240        
        
        
        Cx[0,-3]=-(19/16)*e**2*cosine(iota)*sine(iota)*sine(3*beta)+65/96*e**4*cosine(iota)*sine(iota)*sine(3*beta)+(19*e**6*cosine(iota)*sine(iota)*sine(3*beta))/6144
        Cx[1,-3]=6*e*cosine(iota)*sine(iota)*sine(3*beta)-11*e**3*cosine(iota)*sine(iota)*sine(3*beta)+133/24*e**5*cosine(iota)*sine(iota)*sine(3*beta)
        Cx[2,-3]=-(9/2)*cosine(iota)*sine(iota)*sine(3*beta)+27*e**2*cosine(iota)*sine(iota)*sine(3*beta)-5319/128*e**4*cosine(iota)*sine(iota)*sine(3*beta)+3141/128*e**6*cosine(iota)*sine(iota)*sine(3*beta)
        Cx[3,-3]=-16*e*cosine(iota)*sine(iota)*sine(3*beta)+76*e**3*cosine(iota)*sine(iota)*sine(3*beta)-352/3*e**5*cosine(iota)*sine(iota)*sine(3*beta)
        Cx[4,-3]=-(625/16)*e**2*cosine(iota)*sine(iota)*sine(3*beta)+5625/32*e**4*cosine(iota)*sine(iota)*sine(3*beta)-(1746875*e**6*cosine(iota)*sine(iota)*sine(3*beta))/6144
        Cx[5,-3]=-81*e**3*cosine(iota)*sine(iota)*sine(3*beta)+729/2*e**5*cosine(iota)*sine(iota)*sine(3*beta)
        Cx[6,-3]=-(117649/768)*e**4*cosine(iota)*sine(iota)*sine(3*beta)+(2705927*e**6*cosine(iota)*sine(iota)*sine(3*beta))/3840
        Cx[7,-3]=-(4096/15)*e**5*cosine(iota)*sine(iota)*sine(3*beta)
        Cx[8,-3]=-((4782969*e**6*cosine(iota)*sine(iota)*sine(3*beta))/10240)        
        

        return(Cx[a-1,b])
    
    
    def scross(self, a, b):
        sine = self.sine_
        cosine = self.cosine_
        
        e = self.e
        beta = self.beta
        iota = self.iota
        Sx = np.zeros((10,5))
        
        Sx[0,1]=-(1/16)*e**2*cosine(iota)*cosine(beta)*sine(iota)-5/96*e**4*cosine(iota)*cosine(beta)*sine(iota)-(187*e**6*cosine(iota)*cosine(beta)*sine(iota))/6144
        Sx[1,1]=-(1/24)*e**5*cosine(iota)*cosine(beta)*sine(iota)
        Sx[2,1]=9/256*e**4*cosine(iota)*cosine(beta)*sine(iota)-(63*e**6*cosine(iota)*cosine(beta)*sine(iota))/1280
        Sx[3,1]=1/15*e**5*cosine(iota)*cosine(beta)*sine(iota)
        Sx[4,1]=(625*e**6*cosine(iota)*cosine(beta)*sine(iota))/6144
        
        
        Sx[0,-1]=-(3/2)*cosine(iota)*cosine(beta)*sine(iota)+e**2*cosine(iota)*cosine(beta)*sine(iota)-5/128*e**4*cosine(iota)*cosine(beta)*sine(iota)+13/384*e**6*cosine(iota)*cosine(beta)*sine(iota)
        Sx[1,-1]=-2*e*cosine(iota)*cosine(beta)*sine(iota)+2*e**3*cosine(iota)*cosine(beta)*sine(iota)-11/24*e**5*cosine(iota)*cosine(beta)*sine(iota)
        Sx[2,-1]=-(45/16)*e**2*cosine(iota)*cosine(beta)*sine(iota)+117/32*e**4*cosine(iota)*cosine(beta)*sine(iota)-(2961*e**6*cosine(iota)*cosine(beta)*sine(iota))/2048
        Sx[3,-1]=-4*e**3*cosine(iota)*cosine(beta)*sine(iota)+19/3*e**5*cosine(iota)*cosine(beta)*sine(iota)
        Sx[4,-1]=-(4375/768)*e**4*cosine(iota)*cosine(beta)*sine(iota)+8125/768*e**6*cosine(iota)*cosine(beta)*sine(iota)
        Sx[5,-1]=-(81/10)*e**5*cosine(iota)*cosine(beta)*sine(iota)
        Sx[6,-1]=-((117649*e**6*cosine(iota)*cosine(beta)*sine(iota))/10240)
        
        
        Sx[0,3]=25/256*e**4*cosine(iota)*cosine(3*beta)*sine(iota)+(179*e**6*cosine(iota)*cosine(3*beta)*sine(iota))/3840
        Sx[1,3]=13/120*e**5*cosine(iota)*cosine(3*beta)*sine(iota)
        Sx[2,3]=(1233*e**6*cosine(iota)*cosine(3*beta)*sine(iota))/10240  
        
        
        Sx[0,-3]=19/16*e**2*cosine(iota)*cosine(3*beta)*sine(iota)-65/96*e**4*cosine(iota)*cosine(3*beta)*sine(iota)-(19*e**6*cosine(iota)*cosine(3*beta)*sine(iota))/6144
        Sx[1,-3]=-6*e*cosine(iota)*cosine(3*beta)*sine(iota)+11*e**3*cosine(iota)*cosine(3*beta)*sine(iota)-133/24*e**5*cosine(iota)*cosine(3*beta)*sine(iota)
        Sx[2,-3]=9/2*cosine(iota)*cosine(3*beta)*sine(iota)-27*e**2*cosine(iota)*cosine(3*beta)*sine(iota)+5319/128*e**4*cosine(iota)*cosine(3*beta)*sine(iota)-3141/128*e**6*cosine(iota)*cosine(3*beta)*sine(iota)
        Sx[3,-3]=16*e*cosine(iota)*cosine(3*beta)*sine(iota)-76*e**3*cosine(iota)*cosine(3*beta)*sine(iota)+352/3*e**5*cosine(iota)*cosine(3*beta)*sine(iota)
        Sx[4,-3]=625/16*e**2*cosine(iota)*cosine(3*beta)*sine(iota)-5625/32*e**4*cosine(iota)*cosine(3*beta)*sine(iota)+(1746875*e**6*cosine(iota)*cosine(3*beta)*sine(iota))/6144
        Sx[5,-3]=81*e**3*cosine(iota)*cosine(3*beta)*sine(iota)-729/2*e**5*cosine(iota)*cosine(3*beta)*sine(iota)
        Sx[6,-3]=117649/768*e**4*cosine(iota)*cosine(3*beta)*sine(iota)-(2705927*e**6*cosine(iota)*cosine(3*beta)*sine(iota))/3840
        Sx[7,-3]=4096/15*e**5*cosine(iota)*cosine(3*beta)*sine(iota)
        Sx[8,-3]=(4782969*e**6*cosine(iota)*cosine(3*beta)*sine(iota))/10240        


        return(Sx[a-1,b])
    
    
    
    
    
    
    
    
    
    
    
    