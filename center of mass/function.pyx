### Computes the waveform for moving sources based on GWSurrogate  ###

import numpy as np

from numpy import conj
from numpy import pi

from cython import cmod

from libc.math cimport pow, sqrt, cos, sin, exp
from libc.math cimport M_PI

cdef extern from "complex.h":
    float complex I
    double complex exp(double complex)


# Define constants: speed of light (in m/s), gravitational constant (in m^3/(kg*s^2)), solar mass (in kg) and megaparsec (in m)
c = 2.99e8
G = 6.67e-11
M_sun = 1.99e30
Mpc = 3.09e22


# Define the min and max funstions, and the factorial and binomial coefficient as C-code
cdef int minf(int a, int b):
    if a < b:
        return a
    else:
        return b

cdef int maxf(int a, int b):
    if a > b:
        return a
    else:
        return b

cdef float fact(int n):
    cdef int i
    cdef float f

    if n < 0:
        return 0
    else:
        i, f = 1, 1.0
        while i <= n:
            f = i*f
            i += 1

        return f

cdef float bino(int n, int k):
    if k < 0 or k > n:
        return 0
    else:
        return fact(n)/(fact(k)*fact(n-k))


# Define the spherical harmonics as C-code
cdef double complex Y(int l, int m, float the, float phi):
    cdef double complex hlp
    cdef int k

    hlp = 0

    mink = maxf(0,m+2)
    maxk = minf(l+2,l+m) + 1
    for k in range(mink,maxk):
        hlp += pow(-1,l-k+m+2)*bino(l+2,k)*bino(l-2,k-m-2)*pow(cos(the/2),2*k-m-2)*pow(sin(the/2),2*l-2*k+m+2)*exp(I*m*phi)

    return sqrt(fact(l+m)*fact(l-m)*(2*l+1)/(4*M_PI*fact(l-2)*fact(l+2)))*hlp


# Create the array of spherical harmonics
cdef double complex y[6][11]

cdef createy(float the, float phi):
    cdef int l, m

    for l in range(2,6):
        for m in range(2*l+1):
            mhlp = cmod(m,l+1)-l*(m/(l+1))

            if mhlp < 0:
                mhlp2 = 11 + mhlp
            else:
                mhlp2 = mhlp

            y[l][mhlp2] = Y(l,mhlp,the,phi)



# Pass the surrogate model, the total mass (in M_sun), distance of the source (in Mpc), the sky position (in degree), the velocity magnitude (in km/s), the direction of the velocity (in degree) and the time steps (in s)
# Define the waveform
def WaVe(surrogate_, M_, q_, dis_, the_, phi_, v_mag_, v_the_, v_phi_, times_):

    sur = surrogate_
    M = M_*M_sun
    q = q_
    dis = dis_*Mpc
    the = pi*float(the_)/180
    phi = pi*float(phi_)/180
    v_mag = 1000*float(v_mag_)/c
    v_the = pi*float(v_the_)/180
    v_phi = pi*float(v_phi_)/180
    times = times_

    times = np.array(times)

    # Set the spins of the BHs        
    chiA = [0,0,0]
    chiB = [0,0,0]

    # Compute velocity vector
    v = np.array([sin(v_the)*cos(v_phi),sin(v_the)*sin(v_phi),cos(v_the)])

    # Compute the Doppler shift and the amplitude of the wave
    Dopp = (1+v_mag*(sin(the)*cos(phi)*v[0]+sin(the)*sin(phi)*v[1]+cos(the)*v[2]))/sqrt(1-v_mag**2)
    amp = G*M/(dis*c**2)

    # Transform the times for evaluation in natural units and set f_low to zero so that it can in principle compute all times
    times = c**3*times/(G*M*Dopp)
    f_low = 0

    # 'Get' the wave form
    t, h, dyn = sur(q, chiA, chiB, times=times, f_low=f_low)

    # Compute and write out the spherical harmonics
    createy(the,phi)

    # Calling the function for Coeff
    #Coeff_ = Coeff.Fn2(v_mag, v)
    #c0, cp, cm = Coeff_.Coefficient()
    #c0, cp, cm = Coefficient(v_mag, v)

    # Compute the beamed modes
    # Create an array for the beamed modes
    hb = {}
    for mi in h:
        hb[mi] = np.zeros(len(t), dtype=complex)
        hb[(mi[0],-mi[1])] = np.zeros(len(t), dtype=complex)


    # Compute the beamed modes over time
    for mi in hb:
        # Compute the contributions of each mode
        hlp0 = np.zeros(len(t), dtype=complex)
        for l0 in range(max(2,abs(mi[1])),6):
            try:
                h[(l0,abs(mi[1]))]
            except KeyError:
                pass
            else:
                if mi[1] < 0:
                    hlp0 = hlp0 + (1/4**l0)*(-1)**l0*conj(h[(l0,abs(mi[1]))])*c0[l0][mi[0]][mi[1]]
                else:
                    hlp0 = hlp0 + (1/4**l0)*h[(l0,mi[1])]*c0[l0][mi[0]][mi[1]]

        hlpp = np.zeros(len(t), dtype=complex)
        for lp in range(max(2,abs(mi[1]+1)),6):
            try:
                h[(lp,abs(mi[1]+1))]
            except KeyError:
                pass
            else:
                if mi[1]+1 < 0:
                    hlpp = hlpp + (1/4**lp)*(-1)**lp*conj(h[(lp,abs(mi[1]+1))])*cp[lp][mi[0]][mi[1]]
                else:
                    hlpp = hlpp + (1/4**lp)*h[(lp,mi[1]+1)]*cp[lp][mi[0]][mi[1]]

        hlpm = np.zeros(len(t), dtype=complex)
        for lm in range(max(2,abs(mi[1]-1)),6):
            try:
                h[(lm,abs(mi[1]-1))]
            except KeyError:
                pass
            else:
                if mi[1]-1 < 0:
                    hlpm = hlpm + (1/4**lm)*(-1)**lm*conj(h[(lm,abs(mi[1]-1))])*cp[lp][mi[0]][mi[1]]
                else:
                    hlpm = hlpm + (1/4**lm)*h[(lm,mi[1]-1)]*cm[lm][mi[0]][mi[1]]

        # Sum up the contributions to get the exited mode
        if mi[1] < 0:
            hb[mi] = (-1)**mi[0]*conj(h[(mi[0],abs(mi[1]))]) + (-1)**abs(mi[1])*(1j/4**(mi[0]+1))*(hlp0+hlpp+hlpm)
        else:
            hb[mi] = h[mi] + (-1)**mi[1]*(1j/4**(mi[0]+1))*(hlp0+hlpp+hlpm)

    # Define arrays for the polarizations
    hp = np.zeros(len(t), dtype=float)
    hc = np.zeros(len(t), dtype=float)

    # Compute polarizations in 'SI' units
    hlp = np.zeros(len(t), dtype=complex)
    for mi in hb:
        hlp = hlp + hb[mi]*y[mi[0]][mi[1]]

    hp, hc = amp*hlp.real, amp*hlp.imag


    # Write out the waveforms
    return (hp,hc)
'''
###############################################################
#Coefficients
###############################################################
### Computes coefficients for the modes excited by motion  ###

# Define a funtion to check if the number is bigger than zero, the factorial and binomial coefficient as C-code
cdef int btz(int a):
    if a > 0:
        return a
    else:
        return 0

cdef float fact(int n):
    cdef int i
    cdef float f

    if n < 0:
        return 0
    else:
        i, f = 1, 1.0
        while i <= n:
            f = i*f
            i += 1

        return f

cdef float bino(int n, int k):
    if k < 0 or k > n:
        return 0
    else:
        return fact(n)/(fact(k)*fact(n-k))


# Define the A, F and G functions as C-code
cdef float A(int l, int m):
    if (l-m)*(l+m+1) > 0 and (l+m)*(l-m+1) > 0:
        return sqrt((l-m)*(l+m+1)) - sqrt((l+m)*(l-m+1))
    elif (l-m)*(l+m+1) > 0:
        return sqrt((l-m)*(l+m+1))
    elif (l+m)*(l-m+1) > 0:
        return - sqrt((l+m)*(l-m+1))
    else:
        return 0

cdef float F(int l1, int m1, int l2, int m2):
        if l1 < abs(m1) or l2 < abs(m2):
            return 0
        else:
            return sqrt(fact(l1+m1)*fact(l1-m1)*(2*l1+1)/(fact(l1-2)*fact(l1+2)))*sqrt(fact(l2+m2)*fact(l2-m2)*(2*l2+1)/(fact(l2-2)*fact(l2+2)))

cdef float G(int l1, int m1, int k1, int a1, int b1, int l2, int m2, int k2, int a2, int b2):
    return bino(l1+2,k1)*bino(l1-2,k1-m1-2)*bino(2*l1-2*k1+m1+2,a1)*bino(2*k1-m1-2,b1)*bino(l2-2,k2)*bino(l2+2,k2+m2+2)*bino(2*l2-2*k2-m2-2,a2)*bino(2*k2+m2+2,b2)


# Define mode coefficient functions as C-code
cdef float complex C0(int l1, int l2, int m, float v_mag, float v0, float v1, float v2):
    cdef float hlp
    cdef int k1, a1, b1, k2, a2, b2, u

    hlp = 0
    for k1 in range(btz(m+2),l1+3):
        for a1 in range(2*l1-2*k1+m+3):
            for b1 in range(2*k1-m-1):
                for k2 in range(btz(m-2),l2-1):
                    for a2 in range(2*l2-2*k2-m-1):
                        for b2 in range(2*k2+m+3):
                            u = l1+l2-a1-a2-b1-b2

                            if u == 0:
                                hlp += pow(-1,a1+a2)*m*G(l1,m,k1,a1,b1,l2,m,k2,a2,b2)

    return 2*v_mag*M_PI*I*v2*F(l1,m,l2,m)*hlp

cdef float complex Cp(int l1, int l2, int m, float v_mag, float v0, float v1, float v2):
    cdef float hlp
    cdef int k1, a1, b1, k2, a2, b2, u

    hlp = 0
    for k1 in range(btz(m+3),l1+3):
        for a1 in range(2*l1-2*k1+m+4):
            for b1 in range(2*k1-m-2):
                for k2 in range(btz(m-2),l2-1):
                    for a2 in range(2*l2-2*k2-m-1):
                        for b2 in range(2*k2+m+3):
                            u = l1+l2-a1-a2-b1-b2

                            if abs(u%2) == 1:
                                hlp += pow(-1,a1+a2)*G(l1,m+1,k1,a1,b1,l2,m,k2,a2,b2)*2*A(l1,m+1)/(u*(u**2-4))
                            elif u == 2:
                                hlp += pow(-1,a1+a2)*G(l1,m+1,k1,a1,b1,l2,m,k2,a2,b2)*(m-1)*M_PI/4
                            elif u == -2:
                                hlp += pow(-1,a1+a2)*G(l1,m+1,k1,a1,b1,l2,m,k2,a2,b2)*(m-1)*M_PI/4

    return (I*v0-v1)*v_mag*F(l1,m+1,l2,m)*hlp

cdef float complex Cm(int l1, int l2, int m, float v_mag, float v0, float v1, float v2):
    cdef float hlp
    cdef int k1, a1, b1, k2, a2, b2, u

    hlp = 0
    for k1 in range(btz(m+1),l1+3):
        for a1 in range(2*l1-2*k1+m+2):
            for b1 in range(2*k1-m):
                for k2 in range(btz(m-2),l2-1):
                    for a2 in range(2*l2-2*k2-m-1):
                        for b2 in range(2*k2+m+3):
                            u = l1+l2-a1-a2-b1-b2

                            if abs(u%2) == 1:
                                hlp += pow(-1,a1+a2)*G(l1,m-1,k1,a1,b1,l2,m,k2,a2,b2)*2*A(l1,m-1)/(u*(u**2-4))
                            elif u == 2:
                                hlp += pow(-1,a1+a2)*G(l1,m-1,k1,a1,b1,l2,m,k2,a2,b2)*(m-1)*M_PI/4
                            elif u == -2:
                                hlp += pow(-1,a1+a2)*G(l1,m-1,k1,a1,b1,l2,m,k2,a2,b2)*(m-1)*M_PI/4

    return (I*v0+v1)*v_mag*F(l1,m-1,l2,m)*hlp


# Create the array of coefficients
cdef double complex coeff[3][6][6][11]

cdef createc(float v_mag, float v0, float v1, float v2):
    cdef int l1, l2, m

    for l1 in range(2,6):
        for l2 in range(2,6):
            for m in range(2*l2+1):
                mhlp = cmod(m,l2+1)-l2*(m/(l2+1))

                if mhlp < 0:
                    mhlp2 = 11 + mhlp
                else:
                    mhlp2 = mhlp

                coeff[0][l1][l2][mhlp2] = C0(l1,l2,mhlp,v_mag,v0,v1,v2)
                coeff[1][l1][l2][mhlp2] = Cp(l1,l2,mhlp,v_mag,v0,v1,v2)
                coeff[2][l1][l2][mhlp2] = Cm(l1,l2,mhlp,v_mag,v0,v1,v2)


# Compute coefficients und save them in an array
def Coefficient(v_mag_, v_):
    v_mag = v_mag_
    v = v_

    # Compute the coefficients
    createc(v_mag,v[0],v[1],v[2])

    return (coeff[0],coeff[1],coeff[2])

'''





















