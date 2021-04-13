### Computes coefficients for the modes excited by motion  ###

import numpy as np

from math import sqrt
from math import factorial
from math import pi

# Define constants: speed of light (in m/s), gravitational constant (in m^3/(kg*s^2)), solar mass (in kg) and megaparsec (in m)
c = 2.99e8
G = 6.67e-11
M_sun = 1.99e30
Mpc = 3.09e22


class Fn2:
    def __init__(self, v_mag_, v_):
        
        self.v_mag_ = v_mag_
        self.v_ = v_
    
    # Define delta, binomial coeffcient, A, F and G functions
    def delta_(self,a,b):
        if a == b:
            return(1)
        else:
            return(0)

    def binomial_(self,n,k):
        if n < 0 or k < 0 or n-k < 0:
            return 0
        else:
            return factorial(n)/(factorial(k)*factorial(n-k))

    def A_(self,l,m):
        if (l-m)*(l+m+1) > 0 and (l+m)*(l-m+1) > 0:
            return sqrt((l-m)*(l+m+1)) - sqrt((l+m)*(l-m+1))
        elif (l-m)*(l+m+1) > 0:
            return sqrt((l-m)*(l+m+1))
        elif (l+m)*(l-m+1) > 0:
            return - sqrt((l+m)*(l-m+1))
        else:
            return 0

    def F_(self,l1,m1,l2,m2):
        if l1 < abs(m1) or l2 < abs(m2):
            return 0
        else:
            return sqrt(factorial(l1+m1)*factorial(l1-m1)*(2*l1+1)/(factorial(l1-2)*factorial(l1+2)))*sqrt(factorial(l2+m2)*factorial(l2-m2)*(2*l2+1)/(factorial(l2-2)*factorial(l2+2)))

    def G_(self,l1,m1,k1,a1,b1,l2,m2,k2,a2,b2):
        binomial = self.binomial_
        return binomial(l1+2,k1)*binomial(l1-2,k1-m1-2)*binomial(2*l1-2*k1+m1+2,a1)*binomial(2*k1-m1-2,b1)*binomial(l2-2,k2)*binomial(l2+2,k2+m2+2)*binomial(2*l2-2*k2-m2-2,a2)*binomial(2*k2+m2+2,b2)


    # Define mode coefficient functions
    def C0_(self,l1,l2,m):
        v_mag = self.v_mag_
        v = self.v_
        delta = self.delta_
        A = self.A_
        F = self.F_
        G = self.G_
        
        hlp = 0

        for k1 in range(max(0,m+2),l1+3):
            for a1 in range(2*l1-2*k1+m+3):
                for b1 in range(2*k1-m-1):
                    for k2 in range(max(0,m-2),l2-1):
                        for a2 in range(2*l2-2*k2-m-1):
                            for b2 in range(2*k2+m+3):
                                u = l1+l2-a1-a2-b1-b2
                                if delta(u,0):
                                    hlp = hlp + (-1)**(a1+a2)*m*G(l1,m,k1,a1,b1,l2,m,k2,a2,b2)

        return 2*v_mag*pi*1j*v[2]*F(l1,m,l2,m)*hlp

    def Cp_(self,l1,l2,m):
        v_mag = self.v_mag_
        v = self.v_
        delta = self.delta_
        A = self.A_
        F = self.F_
        G = self.G_
        
        hlp = 0

        for k1 in range(max(0,m+3),l1+3):
            for a1 in range(2*l1-2*k1+m+4):
                for b1 in range(2*k1-m-2):
                    for k2 in range(max(0,m-2),l2-1):
                        for a2 in range(2*l2-2*k2-m-1):
                            for b2 in range(2*k2+m+3):
                                u = l1+l2-a1-a2-b1-b2

                                if delta(abs(u%2),1):
                                    hlp = hlp + (-1)**(a1+a2)*G(l1,m+1,k1,a1,b1,l2,m,k2,a2,b2)*2*A(l1,m+1)/(u*(u**2-4))
                                elif delta(u,2):
                                    hlp = hlp + (-1)**(a1+a2)*G(l1,m+1,k1,a1,b1,l2,m,k2,a2,b2)*(m-1)*pi/4
                                elif delta(u,-2):
                                    hlp = hlp - (-1)**(a1+a2)*G(l1,m+1,k1,a1,b1,l2,m,k2,a2,b2)*(m-1)*pi/4

        return (1j*v[0]-v[1])*v_mag*F(l1,m+1,l2,m)*hlp

    def Cm_(self,l1,l2,m):
        v_mag = self.v_mag_
        v = self.v_
        delta = self.delta_
        A = self.A_
        F = self.F_
        G = self.G_
        
        hlp = 0

        for k1 in range(max(0,m+1),l1+3):
            for a1 in range(2*l1-2*k1+m+2):
                for b1 in range(2*k1-m):
                    for k2 in range(max(0,m-2),l2-1):
                        for a2 in range(2*l2-2*k2-m-1):
                            for b2 in range(2*k2+m+3):
                                u = l1+l2-a1-a2-b1-b2

                                if delta(abs(u%2),1):
                                    hlp = hlp + (-1)**(a1+a2)*G(l1,m-1,k1,a1,b1,l2,m,k2,a2,b2)*2*A(l1,m-1)/(u*(u**2-4))
                                elif delta(u,2):
                                    hlp = hlp - (-1)**(a1+a2)*G(l1,m-1,k1,a1,b1,l2,m,k2,a2,b2)*(m-1)*pi/4
                                elif delta(u,-2):
                                    hlp = hlp + (-1)**(a1+a2)*G(l1,m-1,k1,a1,b1,l2,m,k2,a2,b2)*(m-1)*pi/4

        return (1j*v[0]+v[1])*v_mag*F(l1,m-1,l2,m)*hlp
    
    
    def Coefficient(self):
        C0 = self.C0_
        Cp = self.Cp_
        Cm = self.Cm_
        
        # Set maximal 'l' available
        l = 5
        
        # Define arrays to save the coefficients
        coeff = np.zeros((3,l+1,l+1,2*(l+1)), dtype=complex)

        # Compute and save the coefficients
        for l1 in range(2,l+1):
            for l2 in range(2,l+1):
                for m in range(-l2,l2+1):
                    coeff[0][l1][l2][m] = C0(l1,l2,m)
                    coeff[1][l1][l2][m] = Cp(l1,l2,m)
                    coeff[2][l1][l2][m] = Cm(l1,l2,m)

                    
        return(coeff[0] ,coeff[1] ,coeff[2])

