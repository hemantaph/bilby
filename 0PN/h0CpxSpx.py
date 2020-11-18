import numpy as np

C = 299792458.0

class Fn:
    def __init__(self, e, beta, iota):
        self.e = e
        self.beta = beta
        self.iota = iota
    
    #take out the first row
    def row(self,matrix):
        matrix2 = np.zeros((np.shape(matrix)[0]-1,np.shape(matrix)[1]))
        aa = np.shape(matrix2)[0]
        bb = np.shape(matrix2)[1]

        for i in range(aa):
            for j in range (bb):

                matrix2[i][j] = matrix[i+1][j]

        return(matrix2) 
    

    def cplus(self):
        et = self.e
        beta = self.beta
        iota = self.iota
        Cp = np.zeros((11,3))

        Cp[1,-2] = et**3*( - (13/16)*np.cos(2*beta) - 13/16*np.cos(iota)**2*np.cos(2*beta)) + et**5*( - (5/384)*np.cos(2*beta) - 5/384*np.cos(iota)**2*np.cos(2*beta)) + et*(3/2*np.cos(2*beta) + 3/2*np.cos(iota)**2*np.cos(2*beta))

        Cp[2,-2] = - 2*np.cos(2*beta) - 2*np.cos(iota)**2*np.cos(2*beta) + et**4*( - (23/8)*np.cos(2*beta) - 23/8*np.cos(iota)**2*np.cos(2*beta)) + et**6*(65/144*np.cos(2*beta) + 65/144*np.cos(iota)**2*np.cos(2*beta)) + et**2*(5*np.cos(2*beta) + 5*np.cos(iota)**2*np.cos(2*beta))

        Cp[3,-2] = et**5*( - (963/128)*np.cos(2*beta) - 963/128*np.cos(iota)**2*np.cos(2*beta)) + et*( - (9/2)*np.cos(2*beta) - 9/2*np.cos(iota)**2*np.cos(2*beta)) + et**3*(171/16*np.cos(2*beta) + 171/16*np.cos(iota)**2*np.cos(2*beta))

        Cp[4,-2] = et**6*( - (101/6)*np.cos(2*beta) - 101/6*np.cos(iota)**2*np.cos(2*beta)) + et**2*( - 8*np.cos(2*beta) - 8*np.cos(iota)**2*np.cos(2*beta)) + et**4*(20*np.cos(2*beta) + 20*np.cos(iota)**2*np.cos(2*beta))

        Cp[5,-2] = et**3*( - (625/48)*np.cos(2*beta) - 625/48*np.cos(iota)**2*np.cos(2*beta)) + et**5*(26875/768*np.cos(2*beta) + 26875/768*np.cos(iota)**2*np.cos(2*beta))

        Cp[6,-2] = et**4*( - (81/4)*np.cos(2*beta) - 81/4*np.cos(iota)**2*np.cos(2*beta)) + et**6*(2349/40*np.cos(2*beta) + 2349/40*np.cos(iota)**2*np.cos(2*beta))

        Cp[7,-2] = et**5*( - ((117649*np.cos(2*beta))/3840) - (117649*np.cos(iota)**2*np.cos(2*beta))/3840)

        Cp[8,-2] = et**6*( - (2048/45)*np.cos(2*beta) - 2048/45*np.cos(iota)**2*np.cos(2*beta))

        Cp[9,-2] = 0

        Cp[10,-2] = 0


        Cp[1,0] = et*np.sin(iota)**2 - 1/8*et**3*np.sin(iota)**2 + 1/192*et**5*np.sin(iota)**2

        Cp[2,0] = et**2*np.sin(iota)**2 - 1/3*et**4*np.sin(iota)**2 + 1/24*et**6*np.sin(iota)**2

        Cp[3,0] = 9/8*et**3*np.sin(iota)**2 - 81/128*et**5*np.sin(iota)**2

        Cp[4,0] = 4/3*et**4*np.sin(iota)**2 - 16/15*et**6*np.sin(iota)**2

        Cp[5,0] = 625/384*et**5*np.sin(iota)**2

        Cp[6,0] = 81/40*et**6*np.sin(iota)**2

        Cp[7,0] = 0

        Cp[8,0] = 0

        Cp[9,0] = 0

        Cp[10,0] = 0


        Cp[1,2] = et**5*(47/768*np.cos(2*beta) + 47/768*np.cos(iota)**2*np.cos(2*beta)) + et**3*(7/48*np.cos(2*beta) + 7/48*np.cos(iota)**2*np.cos(2*beta))

        Cp[2,2] = et**6*(11/240*np.cos(2*beta) + 11/240*np.cos(iota)**2*np.cos(2*beta)) + et**4*(1/8*np.cos(2*beta) + 1/8*np.cos(iota)**2*np.cos(2*beta))

        Cp[3,2] = et**5*((153*np.cos(2*beta))/1280 + (153*np.cos(iota)**2*np.cos(2*beta))/1280)

        Cp[4,2] = et**6*(11/90*np.cos(2*beta) + 11/90*np.cos(iota)**2*np.cos(2*beta))

        Cp[5,2] = 0

        Cp[6,2] = 0

        Cp[7,2] = 0

        Cp[8,2] = 0

        Cp[9,2] = 0

        Cp[10,2] = 0

        return(self.row(Cp))
    
    
    def splus(self):
        et = self.e
        beta = self.beta
        iota = self.iota
        Sp = np.zeros((11,3))

        Sp[1,-2] = et**3*( - (13/16)*np.sin(2*beta) - 13/16*np.cos(iota)**2*np.sin(2*beta)) + et**5*( - (5/384)*np.sin(2*beta) - 5/384*np.cos(iota)**2*np.sin(2*beta)) + et*(3/2*np.sin(2*beta) + 3/2*np.cos(iota)**2*np.sin(2*beta))

        Sp[2,-2] = - 2*np.sin(2*beta) - 2*np.cos(iota)**2*np.sin(2*beta) + et**4*( - (23/8)*np.sin(2*beta) - 23/8*np.cos(iota)**2*np.sin(2*beta)) + et**6*(65/144*np.sin(2*beta) + 65/144*np.cos(iota)**2*np.sin(2*beta)) + et**2*(5*np.sin(2*beta) + 5*np.cos(iota)**2*np.sin(2*beta))

        Sp[3,-2] = et**5*( - (963/128)*np.sin(2*beta) - 963/128*np.cos(iota)**2*np.sin(2*beta)) + et*( - (9/2)*np.sin(2*beta) - 9/2*np.cos(iota)**2*np.sin(2*beta)) + et**3*(171/16*np.sin(2*beta) + 171/16*np.cos(iota)**2*np.sin(2*beta))

        Sp[4,-2] = et**6*( - (101/6)*np.sin(2*beta) - 101/6*np.cos(iota)**2*np.sin(2*beta)) + et**2*( - 8*np.sin(2*beta) - 8*np.cos(iota)**2*np.sin(2*beta)) + et**4*(20*np.sin(2*beta) + 20*np.cos(iota)**2*np.sin(2*beta))

        Sp[5,-2] = et**3*( - (625/48)*np.sin(2*beta) - 625/48*np.cos(iota)**2*np.sin(2*beta)) + et**5*(26875/768*np.sin(2*beta) + 26875/768*np.cos(iota)**2*np.sin(2*beta))

        Sp[6,-2] = et**4*( - (81/4)*np.sin(2*beta) - 81/4*np.cos(iota)**2*np.sin(2*beta)) + et**6*(2349/40*np.sin(2*beta) + 2349/40*np.cos(iota)**2*np.sin(2*beta))

        Sp[7,-2] = et**5*( - ((117649*np.sin(2*beta))/3840) - (117649*np.cos(iota)**2*np.sin(2*beta))/3840)

        Sp[8,-2] = et**6*( - (2048/45)*np.sin(2*beta) - 2048/45*np.cos(iota)**2*np.sin(2*beta))

        Sp[9,-2] = 0

        Sp[10,-2] = 0


        Sp[1,0] = 0

        Sp[2,0] = 0

        Sp[3,0] = 0

        Sp[4,0] = 0

        Sp[5,0] = 0

        Sp[6,0] = 0

        Sp[7,0] = 0

        Sp[8,0] = 0

        Sp[9,0] = 0

        Sp[10,0] = 0


        Sp[1,2] = et**3*( - (7/48)*np.sin(2*beta) - 7/48*np.cos(iota)**2*np.sin(2*beta)) + et**5*( - (47/768)*np.sin(2*beta) - 47/768*np.cos(iota)**2*np.sin(2*beta))

        Sp[2,2] = et**4*( - (1/8)*np.sin(2*beta) - 1/8*np.cos(iota)**2*np.sin(2*beta)) + et**6*( - (11/240)*np.sin(2*beta) - 11/240*np.cos(iota)**2*np.sin(2*beta))

        Sp[3,2] = et**5*( - ((153*np.sin(2*beta))/1280) - (153*np.cos(iota)**2*np.sin(2*beta))/1280)

        Sp[4,2] = et**6*( - (11/90)*np.sin(2*beta) - 11/90*np.cos(iota)**2*np.sin(2*beta))

        Sp[5,2] = 0

        Sp[6,2] = 0

        Sp[7,2] = 0

        Sp[8,2] = 0

        Sp[9,2] = 0

        Sp[10,2] = 0

        return(self.row(Sp))
    
    
    def ccross(self):
        et = self.e
        beta = self.beta
        iota = self.iota
        Cx = np.zeros((11,3))

        Cx[1,-2] = - 3*et*np.cos(iota)*np.sin(2*beta) + 13/8*et**3*np.cos(iota)*np.sin(2*beta) + 5/192*et**5*np.cos(iota)*np.sin(2*beta)

        Cx[2,-2] = 4*np.cos(iota)*np.sin(2*beta) - 10*et**2*np.cos(iota)*np.sin(2*beta) + 23/4*et**4*np.cos(iota)*np.sin(2*beta) - 65/72*et**6*np.cos(iota)*np.sin(2*beta)

        Cx[3,-2] = 9*et*np.cos(iota)*np.sin(2*beta) - 171/8*et**3*np.cos(iota)*np.sin(2*beta) + 963/64*et**5*np.cos(iota)*np.sin(2*beta)

        Cx[4,-2] = 16*et**2*np.cos(iota)*np.sin(2*beta) - 40*et**4*np.cos(iota)*np.sin(2*beta) + 101/3*et**6*np.cos(iota)*np.sin(2*beta)

        Cx[5,-2] = 625/24*et**3*np.cos(iota)*np.sin(2*beta) - 26875/384*et**5*np.cos(iota)*np.sin(2*beta)

        Cx[6,-2] = 81/2*et**4*np.cos(iota)*np.sin(2*beta) - 2349/20*et**6*np.cos(iota)*np.sin(2*beta)

        Cx[7,-2] = (117649*et**5*np.cos(iota)*np.sin(2*beta))/1920

        Cx[8,-2] = 4096/45*et**6*np.cos(iota)*np.sin(2*beta)

        Cx[9,-2] = 0

        Cx[10,-2] = 0


        Cx[1,0] = 0

        Cx[2,0] = 0

        Cx[3,0] = 0

        Cx[4,0] = 0

        Cx[5,0] = 0

        Cx[6,0] = 0

        Cx[7,0] = 0

        Cx[8,0] = 0

        Cx[9,0] = 0

        Cx[10,0] = 0


        Cx[1,2] = - (7/24)*et**3*np.cos(iota)*np.sin(2*beta) - 47/384*et**5*np.cos(iota)*np.sin(2*beta)

        Cx[2,2] = - (1/4)*et**4*np.cos(iota)*np.sin(2*beta) - 11/120*et**6*np.cos(iota)*np.sin(2*beta)

        Cx[3,2] = - (153/640)*et**5*np.cos(iota)*np.sin(2*beta)

        Cx[4,2] = - (11/45)*et**6*np.cos(iota)*np.sin(2*beta)

        Cx[5,2] = 0

        Cx[6,2] = 0

        Cx[7,2] = 0

        Cx[8,2] = 0

        Cx[9,2] = 0

        Cx[10,2] = 0

        return(self.row(Cx))
    
    
    def scross(self):
        et = self.e
        beta = self.beta
        iota = self.iota
        Sx = np.zeros((11,3))
        
        Sx[1,-2] = 3*et*np.cos(iota)*np.cos(2*beta) - 13/8*et**3*np.cos(iota)*np.cos(2*beta) - 5/192*et**5*np.cos(iota)*np.cos(2*beta)

        Sx[2,-2] = - 4*np.cos(iota)*np.cos(2*beta) + 10*et**2*np.cos(iota)*np.cos(2*beta) - 23/4*et**4*np.cos(iota)*np.cos(2*beta) + 65/72*et**6*np.cos(iota)*np.cos(2*beta)

        Sx[3,-2] = - 9*et*np.cos(iota)*np.cos(2*beta) + 171/8*et**3*np.cos(iota)*np.cos(2*beta) - 963/64*et**5*np.cos(iota)*np.cos(2*beta)

        Sx[4,-2] = - 16*et**2*np.cos(iota)*np.cos(2*beta) + 40*et**4*np.cos(iota)*np.cos(2*beta) - 101/3*et**6*np.cos(iota)*np.cos(2*beta)

        Sx[5,-2] = - (625/24)*et**3*np.cos(iota)*np.cos(2*beta) + 26875/384*et**5*np.cos(iota)*np.cos(2*beta)

        Sx[6,-2] = - (81/2)*et**4*np.cos(iota)*np.cos(2*beta) + 2349/20*et**6*np.cos(iota)*np.cos(2*beta)

        Sx[7,-2] = - ((117649*et**5*np.cos(iota)*np.cos(2*beta))/1920)

        Sx[8,-2] = - (4096/45)*et**6*np.cos(iota)*np.cos(2*beta)

        Sx[9,-2] = 0

        Sx[10,-2] = 0


        Sx[1,0] = 0

        Sx[2,0] = 0

        Sx[3,0] = 0

        Sx[4,0] = 0

        Sx[5,0] = 0

        Sx[6,0] = 0

        Sx[7,0] = 0

        Sx[8,0] = 0

        Sx[9,0] = 0

        Sx[10,0] = 0


        Sx[1,2] = - (7/24)*et**3*np.cos(iota)*np.cos(2*beta) - 47/384*et**5*np.cos(iota)*np.cos(2*beta)

        Sx[2,2] = - (1/4)*et**4*np.cos(iota)*np.cos(2*beta) - 11/120*et**6*np.cos(iota)*np.cos(2*beta)

        Sx[3,2] = - (153/640)*et**5*np.cos(iota)*np.cos(2*beta)

        Sx[4,2] = - (11/45)*et**6*np.cos(iota)*np.cos(2*beta)

        Sx[5,2] = 0

        Sx[6,2] = 0

        Sx[7,2] = 0

        Sx[8,2] = 0

        Sx[9,2] = 0

        Sx[10,2] = 0

        return(self.row(Sx))
    
    
    
    
    
    
    
    
    
    
    
    