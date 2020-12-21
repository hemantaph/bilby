import numpy as np
import sympy as sp

C = 299792458.0

class Fn:
    def __init__(self, e, beta, iota):      
        self.e = e
        self.beta = beta
        self.iota = iota
    

    def cplus(self, a, b):
        e = self.e
        beta = self.beta
        iota = self.iota
        Cp = np.zeros((10,5))

        Cp[0,1]=e**6 * ((109*sp.cos(beta)*sp.sin(iota))/12288+(265*sp.cos(iota)**2*sp.cos(beta)*sp.sin(iota))/12288)+e**4 * (1/64*sp.cos(beta)*sp.sin(iota)+7/192*sp.cos(iota)**2*sp.cos(beta)*sp.sin(iota))+e**2 * (-(9/32)*sp.cos(beta)*sp.sin(iota)+11/32*sp.cos(iota)**2*sp.cos(beta)*sp.sin(iota))
        Cp[1,1]=e**5 * (1/16*sp.cos(beta)*sp.sin(iota)-1/48*sp.cos(iota)**2*sp.cos(beta)*sp.sin(iota))+e**3 * (-(1/3)*sp.cos(beta)*sp.sin(iota)+1/3*sp.cos(iota)**2*sp.cos(beta)*sp.sin(iota))
        Cp[2,1]=e**6 * ((369*sp.cos(beta)*sp.sin(iota))/2560-(243*sp.cos(iota)**2*sp.cos(beta)*sp.sin(iota))/2560)+e**4 * (-(207/512)*sp.cos(beta)*sp.sin(iota)+189/512*sp.cos(iota)**2*sp.cos(beta)*sp.sin(iota))
        Cp[3,1]=e**5 * (-(1/2)*sp.cos(beta)*sp.sin(iota)+13/30*sp.cos(iota)**2*sp.cos(beta)*sp.sin(iota))
        Cp[4,1]=e**6 * (-((23125*sp.cos(beta)*sp.sin(iota))/36864)+(19375*sp.cos(iota)**2*sp.cos(beta)*sp.sin(iota))/36864)


        Cp[0,-1]=-(5/4)*sp.cos(beta)*sp.sin(iota)-1/4*sp.cos(iota)**2*sp.cos(beta)*sp.sin(iota)+e**2 * (3/2*sp.cos(beta)*sp.sin(iota)-1/2*sp.cos(iota)**2*sp.cos(beta)*sp.sin(iota))+e**6 * ((41*sp.cos(beta)*sp.sin(iota))/2304+(37*sp.cos(iota)**2*sp.cos(beta)*sp.sin(iota))/2304)+e**4 * (-(51/256)*sp.cos(beta)*sp.sin(iota)+41/256*sp.cos(iota)**2*sp.cos(beta)*sp.sin(iota))
        Cp[1,-1]=e**3 * (4*sp.cos(beta)*sp.sin(iota)-2*sp.cos(iota)**2*sp.cos(beta)*sp.sin(iota))+e**5 * (-(19/16)*sp.cos(beta)*sp.sin(iota)+35/48*sp.cos(iota)**2*sp.cos(beta)*sp.sin(iota))+e * (-3*sp.cos(beta)*sp.sin(iota)+sp.cos(iota)**2*sp.cos(beta)*sp.sin(iota))
        Cp[2,-1]=e**4 * (531/64*sp.cos(beta)*sp.sin(iota)-297/64*sp.cos(iota)**2*sp.cos(beta)*sp.sin(iota))+e**6 * (-((15399*sp.cos(beta)*sp.sin(iota))/4096)+(9477*sp.cos(iota)**2*sp.cos(beta)*sp.sin(iota))/4096)+e**2 * (-(171/32)*sp.cos(beta)*sp.sin(iota)+81/32*sp.cos(iota)**2*sp.cos(beta)*sp.sin(iota))
        Cp[3,-1]=e**5 * (31/2*sp.cos(beta)*sp.sin(iota)-55/6*sp.cos(iota)**2*sp.cos(beta)*sp.sin(iota))+e**3 * (-(26/3)*sp.cos(beta)*sp.sin(iota)+14/3*sp.cos(iota)**2*sp.cos(beta)*sp.sin(iota))
        Cp[4,-1]=e**6 * ((41875*sp.cos(beta)*sp.sin(iota))/1536-(25625*sp.cos(iota)**2*sp.cos(beta)*sp.sin(iota))/1536)+e**4 * (-(6875/512)*sp.cos(beta)*sp.sin(iota)+(11875*sp.cos(iota)**2*sp.cos(beta)*sp.sin(iota))/1536)
        Cp[5,-1]=e**5 * (-(81/4)*sp.cos(beta)*sp.sin(iota)+243/20*sp.cos(iota)**2*sp.cos(beta)*sp.sin(iota))
        Cp[6,-1]=e**6 * (-((5529503*sp.cos(beta)*sp.sin(iota))/184320)+(3411821*sp.cos(iota)**2*sp.cos(beta)*sp.sin(iota))/184320)
        
        
        Cp[0,3]=e**4 * (-(25/512)*sp.cos(3*beta)*sp.sin(iota)-25/512*sp.cos(iota)**2*sp.cos(3*beta)*sp.sin(iota))+e**6 * (-((179*sp.cos(3*beta)*sp.sin(iota))/7680)-(179*sp.cos(iota)**2*sp.cos(3*beta)*sp.sin(iota))/7680)
        Cp[1,3]=e**5 * (-(13/240)*sp.cos(3*beta)*sp.sin(iota)-13/240*sp.cos(iota)**2*sp.cos(3*beta)*sp.sin(iota))
        Cp[2,3]=e**6 * (-((1233*sp.cos(3*beta)*sp.sin(iota))/20480)-(1233*sp.cos(iota)**2*sp.cos(3*beta)*sp.sin(iota))/20480)
      
    
        Cp[0,-3]=e**4 * (-(65/192)*sp.cos(3*beta)*sp.sin(iota)-65/192*sp.cos(iota)**2*sp.cos(3*beta)*sp.sin(iota))+e**6 * (-((19*sp.cos(3*beta)*sp.sin(iota))/12288)-(19*sp.cos(iota)**2*sp.cos(3*beta)*sp.sin(iota))/12288)+e**2 * (19/32*sp.cos(3*beta)*sp.sin(iota)+19/32*sp.cos(iota)**2*sp.cos(3*beta)*sp.sin(iota))
        Cp[1,-3]=e * (-3*sp.cos(3*beta)*sp.sin(iota)-3*sp.cos(iota)**2*sp.cos(3*beta)*sp.sin(iota))+e**5 * (-(133/48)*sp.cos(3*beta)*sp.sin(iota)-133/48*sp.cos(iota)**2*sp.cos(3*beta)*sp.sin(iota))+e**3 * (11/2*sp.cos(3*beta)*sp.sin(iota)+11/2*sp.cos(iota)**2*sp.cos(3*beta)*sp.sin(iota))
        Cp[2,-3]=9/4*sp.cos(3*beta)*sp.sin(iota)+9/4*sp.cos(iota)**2*sp.cos(3*beta)*sp.sin(iota)+e**2 * (-(27/2)*sp.cos(3*beta)*sp.sin(iota)-27/2*sp.cos(iota)**2*sp.cos(3*beta)*sp.sin(iota))+e**6 * (-(3141/256)*sp.cos(3*beta)*sp.sin(iota)-3141/256*sp.cos(iota)**2*sp.cos(3*beta)*sp.sin(iota))+e**4 * (5319/256*sp.cos(3*beta)*sp.sin(iota)+5319/256*sp.cos(iota)**2*sp.cos(3*beta)*sp.sin(iota))
        Cp[3,-3]=e**3 * (-38*sp.cos(3*beta)*sp.sin(iota)-38*sp.cos(iota)**2*sp.cos(3*beta)*sp.sin(iota))+e * (8*sp.cos(3*beta)*sp.sin(iota)+8*sp.cos(iota)**2*sp.cos(3*beta)*sp.sin(iota))+e**5 * (176/3*sp.cos(3*beta)*sp.sin(iota)+176/3*sp.cos(iota)**2*sp.cos(3*beta)*sp.sin(iota))
        Cp[4,-3]=e**4 * (-(5625/64)*sp.cos(3*beta)*sp.sin(iota)-5625/64*sp.cos(iota)**2*sp.cos(3*beta)*sp.sin(iota))+e**2 * (625/32*sp.cos(3*beta)*sp.sin(iota)+625/32*sp.cos(iota)**2*sp.cos(3*beta)*sp.sin(iota))+e**6 * ((1746875*sp.cos(3*beta)*sp.sin(iota))/12288+(1746875*sp.cos(iota)**2*sp.cos(3*beta)*sp.sin(iota))/12288)
        Cp[5,-3]=e**5 * (-(729/4)*sp.cos(3*beta)*sp.sin(iota)-729/4*sp.cos(iota)**2*sp.cos(3*beta)*sp.sin(iota))+e**3 * (81/2*sp.cos(3*beta)*sp.sin(iota)+81/2*sp.cos(iota)**2*sp.cos(3*beta)*sp.sin(iota))
        Cp[6,-3]=e**6 * (-((2705927*sp.cos(3*beta)*sp.sin(iota))/7680)-(2705927*sp.cos(iota)**2*sp.cos(3*beta)*sp.sin(iota))/7680)+e**4 * ((117649*sp.cos(3*beta)*sp.sin(iota))/1536+(117649*sp.cos(iota)**2*sp.cos(3*beta)*sp.sin(iota))/1536)
        Cp[7,-3]=e**5 * (2048/15*sp.cos(3*beta)*sp.sin(iota)+2048/15*sp.cos(iota)**2*sp.cos(3*beta)*sp.sin(iota))
        Cp[8,-3]=e**6 * ((4782969*sp.cos(3*beta)*sp.sin(iota))/20480+(4782969*sp.cos(iota)**2*sp.cos(3*beta)*sp.sin(iota))/20480)  


        return(Cp[a-1,b])
    
    
    def splus(self, a, b):
        e = self.e
        beta = self.beta
        iota = self.iota
        Sp = np.zeros((10,5))

        Sp[0,-1] =  - (5/4)*sp.sin(iota)*sp.sin(beta) - 1/4*sp.cos(iota)**2*sp.sin(iota)*sp.sin(beta) + e**2 * (3/2*sp.sin(iota)*sp.sin(beta) - 1/2*sp.cos(iota)**2*sp.sin(iota)*sp.sin(beta)) + e**6 * ((41*sp.sin(iota)*sp.sin(beta))/2304 + (37*sp.cos(iota)**2*sp.sin(iota)*sp.sin(beta))/2304) + e**4 * ( - (51/256)*sp.sin(iota)*sp.sin(beta) + 41/256*sp.cos(iota)**2*sp.sin(iota)*sp.sin(beta))

        Sp[1,-1] = e**3 * (4*sp.sin(iota)*sp.sin(beta) - 2*sp.cos(iota)**2*sp.sin(iota)*sp.sin(beta)) + e**5 * ( - (19/16)*sp.sin(iota)*sp.sin(beta) + 35/48*sp.cos(iota)**2*sp.sin(iota)*sp.sin(beta)) + e * ( - 3*sp.sin(iota)*sp.sin(beta) + sp.cos(iota)**2*sp.sin(iota)*sp.sin(beta))

        Sp[2,-1] = e**4 * (531/64*sp.sin(iota)*sp.sin(beta) - 297/64*sp.cos(iota)**2*sp.sin(iota)*sp.sin(beta)) + e**6 * ( - ((15399*sp.sin(iota)*sp.sin(beta))/4096) + (9477*sp.cos(iota)**2*sp.sin(iota)*sp.sin(beta))/4096) + e**2 * ( - (171/32)*sp.sin(iota)*sp.sin(beta) + 81/32*sp.cos(iota)**2*sp.sin(iota)*sp.sin(beta))

        Sp[3,-1] = e**5 * (31/2*sp.sin(iota)*sp.sin(beta) - 55/6*sp.cos(iota)**2*sp.sin(iota)*sp.sin(beta)) + e**3 * ( - (26/3)*sp.sin(iota)*sp.sin(beta) + 14/3*sp.cos(iota)**2*sp.sin(iota)*sp.sin(beta))

        Sp[4,-1] = e**6 * ((41875*sp.sin(iota)*sp.sin(beta))/1536 - (25625*sp.cos(iota)**2*sp.sin(iota)*sp.sin(beta))/1536) + e**4 * ( - (6875/512)*sp.sin(iota)*sp.sin(beta) + (11875*sp.cos(iota)**2*sp.sin(iota)*sp.sin(beta))/1536)

        Sp[5,-1] = e**5 * ( - (81/4)*sp.sin(iota)*sp.sin(beta) + 243/20*sp.cos(iota)**2*sp.sin(iota)*sp.sin(beta))

        Sp[6,-1] = e**6 * ( - ((5529503*sp.sin(iota)*sp.sin(beta))/184320) + (3411821*sp.cos(iota)**2*sp.sin(iota)*sp.sin(beta))/184320)


        Sp[0,1] = e**2 * (9/32*sp.sin(iota)*sp.sin(beta) - 11/32*sp.cos(iota)**2*sp.sin(iota)*sp.sin(beta)) + e**4 * ( - (1/64)*sp.sin(iota)*sp.sin(beta) - 7/192*sp.cos(iota)**2*sp.sin(iota)*sp.sin(beta)) + e**6 * ( - ((109*sp.sin(iota)*sp.sin(beta))/12288) - (265*sp.cos(iota)**2*sp.sin(iota)*sp.sin(beta))/12288)

        Sp[1,1] = e**3 * (1/3*sp.sin(iota)*sp.sin(beta) - 1/3*sp.cos(iota)**2*sp.sin(iota)*sp.sin(beta)) + e**5 * ( - (1/16)*sp.sin(iota)*sp.sin(beta) + 1/48*sp.cos(iota)**2*sp.sin(iota)*sp.sin(beta))

        Sp[2,1] = e**4 * (207/512*sp.sin(iota)*sp.sin(beta) - 189/512*sp.cos(iota)**2*sp.sin(iota)*sp.sin(beta)) + e**6 * ( - ((369*sp.sin(iota)*sp.sin(beta))/2560) + (243*sp.cos(iota)**2*sp.sin(iota)*sp.sin(beta))/2560)

        Sp[3,1] = e**5 * (1/2*sp.sin(iota)*sp.sin(beta) - 13/30*sp.cos(iota)**2*sp.sin(iota)*sp.sin(beta))

        Sp[4,1] = e**6 * ((23125*sp.sin(iota)*sp.sin(beta))/36864 - (19375*sp.cos(iota)**2*sp.sin(iota)*sp.sin(beta))/36864)


        Sp[0,-3] = e**4 * ( - (65/192)*sp.sin(iota)*sp.sin(3*beta) - 65/192*sp.cos(iota)**2*sp.sin(iota)*sp.sin(3*beta)) + e**6 * ( - ((19*sp.sin(iota)*sp.sin(3*beta))/12288) - (19*sp.cos(iota)**2*sp.sin(iota)*sp.sin(3*beta))/12288) + e**2 * (19/32*sp.sin(iota)*sp.sin(3*beta) + 19/32*sp.cos(iota)**2*sp.sin(iota)*sp.sin(3*beta))

        Sp[1,-3] = e * ( - 3*sp.sin(iota)*sp.sin(3*beta) - 3*sp.cos(iota)**2*sp.sin(iota)*sp.sin(3*beta)) + e**5 * ( - (133/48)*sp.sin(iota)*sp.sin(3*beta) - 133/48*sp.cos(iota)**2*sp.sin(iota)*sp.sin(3*beta)) + e**3 * (11/2*sp.sin(iota)*sp.sin(3*beta) + 11/2*sp.cos(iota)**2*sp.sin(iota)*sp.sin(3*beta))

        Sp[2,-3] = 9/4*sp.sin(iota)*sp.sin(3*beta) + 9/4*sp.cos(iota)**2*sp.sin(iota)*sp.sin(3*beta) + e**2 * ( - (27/2)*sp.sin(iota)*sp.sin(3*beta) - 27/2*sp.cos(iota)**2*sp.sin(iota)*sp.sin(3*beta)) + e**6 * ( - (3141/256)*sp.sin(iota)*sp.sin(3*beta) - 3141/256*sp.cos(iota)**2*sp.sin(iota)*sp.sin(3*beta)) + e**4 * (5319/256*sp.sin(iota)*sp.sin(3*beta) + 5319/256*sp.cos(iota)**2*sp.sin(iota)*sp.sin(3*beta))

        Sp[3,-3] = e**3 * ( - 38*sp.sin(iota)*sp.sin(3*beta) - 38*sp.cos(iota)**2*sp.sin(iota)*sp.sin(3*beta)) + e * (8*sp.sin(iota)*sp.sin(3*beta) + 8*sp.cos(iota)**2*sp.sin(iota)*sp.sin(3*beta)) + e**5 * (176/3*sp.sin(iota)*sp.sin(3*beta) + 176/3*sp.cos(iota)**2*sp.sin(iota)*sp.sin(3*beta))

        Sp[4,-3] = e**4 * ( - (5625/64)*sp.sin(iota)*sp.sin(3*beta) - 5625/64*sp.cos(iota)**2*sp.sin(iota)*sp.sin(3*beta)) + e**2 * (625/32*sp.sin(iota)*sp.sin(3*beta) + 625/32*sp.cos(iota)**2*sp.sin(iota)*sp.sin(3*beta)) + e**6 * ((1746875*sp.sin(iota)*sp.sin(3*beta))/12288 + (1746875*sp.cos(iota)**2*sp.sin(iota)*sp.sin(3*beta))/12288)

        Sp[5,-3] = e**5 * ( - (729/4)*sp.sin(iota)*sp.sin(3*beta) - 729/4*sp.cos(iota)**2*sp.sin(iota)*sp.sin(3*beta)) + e**3 * (81/2*sp.sin(iota)*sp.sin(3*beta) + 81/2*sp.cos(iota)**2*sp.sin(iota)*sp.sin(3*beta))

        Sp[6,-3] = e**6 * ( - ((2705927*sp.sin(iota)*sp.sin(3*beta))/7680) - (2705927*sp.cos(iota)**2*sp.sin(iota)*sp.sin(3*beta))/7680) + e**4 * ((117649*sp.sin(iota)*sp.sin(3*beta))/1536 + (117649*sp.cos(iota)**2*sp.sin(iota)*sp.sin(3*beta))/1536)

        Sp[7,-3] = e**5 * (2048/15*sp.sin(iota)*sp.sin(3*beta) + 2048/15*sp.cos(iota)**2*sp.sin(iota)*sp.sin(3*beta))

        Sp[8,-3] = e**6 * ((4782969*sp.sin(iota)*sp.sin(3*beta))/20480 + (4782969*sp.cos(iota)**2*sp.sin(iota)*sp.sin(3*beta))/20480)


        Sp[0,3] = e**6 * ((179*sp.sin(iota)*sp.sin(3*beta))/7680 + (179*sp.cos(iota)**2*sp.sin(iota)*sp.sin(3*beta))/7680) + e**4 * (25/512*sp.sin(iota)*sp.sin(3*beta) + 25/512*sp.cos(iota)**2*sp.sin(iota)*sp.sin(3*beta))

        Sp[1,3] = e**5 * (13/240*sp.sin(iota)*sp.sin(3*beta) + 13/240*sp.cos(iota)**2*sp.sin(iota)*sp.sin(3*beta))

        Sp[2,3] = e**6 * ((1233*sp.sin(iota)*sp.sin(3*beta))/20480 + (1233*sp.cos(iota)**2*sp.sin(iota)*sp.sin(3*beta))/20480)

        return(Sp[a-1,b])
    
    
    def ccross(self, a, b):
        e = self.e
        beta = self.beta
        iota = self.iota
        Cx = np.zeros((10,5))

        Cx[0,1]=-(1/16)*e**2 *sp.cos(iota)*sp.sin(iota)*sp.sin(beta)-5/96*e**4 *sp.cos(iota)*sp.sin(iota)*sp.sin(beta)-(187*e**6 *sp.cos(iota)*sp.sin(iota)*sp.sin(beta))/6144
        Cx[1,1]=-(1/24)*e**5 *sp.cos(iota)*sp.sin(iota)*sp.sin(beta)
        Cx[2,1]=9/256*e**4 *sp.cos(iota)*sp.sin(iota)*sp.sin(beta)-(63*e**6 *sp.cos(iota)*sp.sin(iota)*sp.sin(beta))/1280
        Cx[3,1]=1/15*e**5 *sp.cos(iota)*sp.sin(iota)*sp.sin(beta)
        Cx[4,1]=(625*e**6 *sp.cos(iota)*sp.sin(iota)*sp.sin(beta))/6144
        
        
        Cx[0,-1]=3/2*sp.cos(iota)*sp.sin(iota)*sp.sin(beta)-e**2 *sp.cos(iota)*sp.sin(iota)*sp.sin(beta)+5/128*e**4 *sp.cos(iota)*sp.sin(iota)*sp.sin(beta)-13/384*e**6 *sp.cos(iota)*sp.sin(iota)*sp.sin(beta)
        Cx[1,-1]=2*e *sp.cos(iota)*sp.sin(iota)*sp.sin(beta)-2*e**3 *sp.cos(iota)*sp.sin(iota)*sp.sin(beta)+11/24*e**5 *sp.cos(iota)*sp.sin(iota)*sp.sin(beta)
        Cx[2,-1]=45/16*e**2 *sp.cos(iota)*sp.sin(iota)*sp.sin(beta)-117/32*e**4 *sp.cos(iota)*sp.sin(iota)*sp.sin(beta)+(2961*e**6 *sp.cos(iota)*sp.sin(iota)*sp.sin(beta))/2048
        Cx[3,-1]=4*e**3 *sp.cos(iota)*sp.sin(iota)*sp.sin(beta)-19/3*e**5 *sp.cos(iota)*sp.sin(iota)*sp.sin(beta)
        Cx[4,-1]=4375/768*e**4 *sp.cos(iota)*sp.sin(iota)*sp.sin(beta)-8125/768*e**6 *sp.cos(iota)*sp.sin(iota)*sp.sin(beta)
        Cx[5,-1]=81/10*e**5 *sp.cos(iota)*sp.sin(iota)*sp.sin(beta)
        Cx[6,-1]=(117649*e**6 *sp.cos(iota)*sp.sin(iota)*sp.sin(beta))/10240 
        
        
        Cx[0,3]=25/256*e**4 *sp.cos(iota)*sp.sin(iota)*sp.sin(3*beta)+(179*e**6 *sp.cos(iota)*sp.sin(iota)*sp.sin(3*beta))/3840
        Cx[1,3]=13/120*e**5 *sp.cos(iota)*sp.sin(iota)*sp.sin(3*beta)
        Cx[2,3]=(1233*e**6 *sp.cos(iota)*sp.sin(iota)*sp.sin(3*beta))/10240        
        
        
        Cx[0,-3]=-(19/16)*e**2 *sp.cos(iota)*sp.sin(iota)*sp.sin(3*beta)+65/96*e**4 *sp.cos(iota)*sp.sin(iota)*sp.sin(3*beta)+(19*e**6 *sp.cos(iota)*sp.sin(iota)*sp.sin(3*beta))/6144
        Cx[1,-3]=6*e *sp.cos(iota)*sp.sin(iota)*sp.sin(3*beta)-11*e**3 *sp.cos(iota)*sp.sin(iota)*sp.sin(3*beta)+133/24*e**5 *sp.cos(iota)*sp.sin(iota)*sp.sin(3*beta)
        Cx[2,-3]=-(9/2)*sp.cos(iota)*sp.sin(iota)*sp.sin(3*beta)+27*e**2 *sp.cos(iota)*sp.sin(iota)*sp.sin(3*beta)-5319/128*e**4 *sp.cos(iota)*sp.sin(iota)*sp.sin(3*beta)+3141/128*e**6 *sp.cos(iota)*sp.sin(iota)*sp.sin(3*beta)
        Cx[3,-3]=-16*e *sp.cos(iota)*sp.sin(iota)*sp.sin(3*beta)+76*e**3 *sp.cos(iota)*sp.sin(iota)*sp.sin(3*beta)-352/3*e**5 *sp.cos(iota)*sp.sin(iota)*sp.sin(3*beta)
        Cx[4,-3]=-(625/16)*e**2 *sp.cos(iota)*sp.sin(iota)*sp.sin(3*beta)+5625/32*e**4 *sp.cos(iota)*sp.sin(iota)*sp.sin(3*beta)-(1746875*e**6 *sp.cos(iota)*sp.sin(iota)*sp.sin(3*beta))/6144
        Cx[5,-3]=-81*e**3 *sp.cos(iota)*sp.sin(iota)*sp.sin(3*beta)+729/2*e**5 *sp.cos(iota)*sp.sin(iota)*sp.sin(3*beta)
        Cx[6,-3]=-(117649/768)*e**4 *sp.cos(iota)*sp.sin(iota)*sp.sin(3*beta)+(2705927*e**6 *sp.cos(iota)*sp.sin(iota)*sp.sin(3*beta))/3840
        Cx[7,-3]=-(4096/15)*e**5 *sp.cos(iota)*sp.sin(iota)*sp.sin(3*beta)
        Cx[8,-3]=-((4782969*e**6 *sp.cos(iota)*sp.sin(iota)*sp.sin(3*beta))/10240)        
        

        return(Cx[a-1,b])
    
    
    def scross(self, a, b):
        e = self.e
        beta = self.beta
        iota = self.iota
        Sx = np.zeros((10,5))
        
        Sx[0,1]=-(1/16)*e**2 *sp.cos(iota)*sp.cos(beta)*sp.sin(iota)-5/96*e**4 *sp.cos(iota)*sp.cos(beta)*sp.sin(iota)-(187*e**6 *sp.cos(iota)*sp.cos(beta)*sp.sin(iota))/6144
        Sx[1,1]=-(1/24)*e**5 *sp.cos(iota)*sp.cos(beta)*sp.sin(iota)
        Sx[2,1]=9/256*e**4 *sp.cos(iota)*sp.cos(beta)*sp.sin(iota)-(63*e**6 *sp.cos(iota)*sp.cos(beta)*sp.sin(iota))/1280
        Sx[3,1]=1/15*e**5 *sp.cos(iota)*sp.cos(beta)*sp.sin(iota)
        Sx[4,1]=(625*e**6 *sp.cos(iota)*sp.cos(beta)*sp.sin(iota))/6144
        
        
        Sx[0,-1]=-(3/2)*sp.cos(iota)*sp.cos(beta)*sp.sin(iota)+e**2 *sp.cos(iota)*sp.cos(beta)*sp.sin(iota)-5/128*e**4 *sp.cos(iota)*sp.cos(beta)*sp.sin(iota)+13/384*e**6 *sp.cos(iota)*sp.cos(beta)*sp.sin(iota)
        Sx[1,-1]=-2*e *sp.cos(iota)*sp.cos(beta)*sp.sin(iota)+2*e**3 *sp.cos(iota)*sp.cos(beta)*sp.sin(iota)-11/24*e**5 *sp.cos(iota)*sp.cos(beta)*sp.sin(iota)
        Sx[2,-1]=-(45/16)*e**2 *sp.cos(iota)*sp.cos(beta)*sp.sin(iota)+117/32*e**4 *sp.cos(iota)*sp.cos(beta)*sp.sin(iota)-(2961*e**6 *sp.cos(iota)*sp.cos(beta)*sp.sin(iota))/2048
        Sx[3,-1]=-4*e**3 *sp.cos(iota)*sp.cos(beta)*sp.sin(iota)+19/3*e**5 *sp.cos(iota)*sp.cos(beta)*sp.sin(iota)
        Sx[4,-1]=-(4375/768)*e**4 *sp.cos(iota)*sp.cos(beta)*sp.sin(iota)+8125/768*e**6 *sp.cos(iota)*sp.cos(beta)*sp.sin(iota)
        Sx[5,-1]=-(81/10)*e**5 *sp.cos(iota)*sp.cos(beta)*sp.sin(iota)
        Sx[6,-1]=-((117649*e**6 *sp.cos(iota)*sp.cos(beta)*sp.sin(iota))/10240)
        
        
        Sx[0,3]=25/256*e**4 *sp.cos(iota)*sp.cos(3*beta)*sp.sin(iota)+(179*e**6 *sp.cos(iota)*sp.cos(3*beta)*sp.sin(iota))/3840
        Sx[1,3]=13/120*e**5 *sp.cos(iota)*sp.cos(3*beta)*sp.sin(iota)
        Sx[2,3]=(1233*e**6 *sp.cos(iota)*sp.cos(3*beta)*sp.sin(iota))/10240  
        
        
        Sx[0,-3]=19/16*e**2 *sp.cos(iota)*sp.cos(3*beta)*sp.sin(iota)-65/96*e**4 *sp.cos(iota)*sp.cos(3*beta)*sp.sin(iota)-(19*e**6 *sp.cos(iota)*sp.cos(3*beta)*sp.sin(iota))/6144
        Sx[1,-3]=-6*e *sp.cos(iota)*sp.cos(3*beta)*sp.sin(iota)+11*e**3 *sp.cos(iota)*sp.cos(3*beta)*sp.sin(iota)-133/24*e**5 *sp.cos(iota)*sp.cos(3*beta)*sp.sin(iota)
        Sx[2,-3]=9/2*sp.cos(iota)*sp.cos(3*beta)*sp.sin(iota)-27*e**2 *sp.cos(iota)*sp.cos(3*beta)*sp.sin(iota)+5319/128*e**4 *sp.cos(iota)*sp.cos(3*beta)*sp.sin(iota)-3141/128*e**6 *sp.cos(iota)*sp.cos(3*beta)*sp.sin(iota)
        Sx[3,-3]=16*e *sp.cos(iota)*sp.cos(3*beta)*sp.sin(iota)-76*e**3 *sp.cos(iota)*sp.cos(3*beta)*sp.sin(iota)+352/3*e**5 *sp.cos(iota)*sp.cos(3*beta)*sp.sin(iota)
        Sx[4,-3]=625/16*e**2 *sp.cos(iota)*sp.cos(3*beta)*sp.sin(iota)-5625/32*e**4 *sp.cos(iota)*sp.cos(3*beta)*sp.sin(iota)+(1746875*e**6 *sp.cos(iota)*sp.cos(3*beta)*sp.sin(iota))/6144
        Sx[5,-3]=81*e**3 *sp.cos(iota)*sp.cos(3*beta)*sp.sin(iota)-729/2*e**5 *sp.cos(iota)*sp.cos(3*beta)*sp.sin(iota)
        Sx[6,-3]=117649/768*e**4 *sp.cos(iota)*sp.cos(3*beta)*sp.sin(iota)-(2705927*e**6 *sp.cos(iota)*sp.cos(3*beta)*sp.sin(iota))/3840
        Sx[7,-3]=4096/15*e**5 *sp.cos(iota)*sp.cos(3*beta)*sp.sin(iota)
        Sx[8,-3]=(4782969*e**6 *sp.cos(iota)*sp.cos(3*beta)*sp.sin(iota))/10240        


        return(Sx[a-1,b])
    
    
    
    
    
    
    
    
    
    
    
    