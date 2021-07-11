### Computes the waveform for moving sources based on GWSurrogate  ###

import gwsurrogate
import Coeff
import numpy as np

from numpy import conj
from numpy import pi

from cython import cmod

from libc.math cimport pow, sqrt, cos, sin, e
from libc.math cimport M_PI

cdef extern from "complex.h":
	double complex I
	double complex cexp(double complex)

import time as time


# Define constants: speed of light (in km/s), gravitational constant (in m^3/(kg*s^2)), solar mass (in kg) and megaparsec (in m)
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
	cdef float f = 1.0

	for i in range(n):
			f = (i+1)*f
			
	return f

cdef float bino(int n, int k):
	if k < 0 or k > n:
		return 0
	else:
		return fact(n)/(fact(k)*fact(n-k))


# Define the spherical harmonics as C-code
cdef double complex Y(int l, int m, float the, float phi):
	cdef double hlp = 0.0
	cdef int k

	mink = maxf(0,m+2)
	maxk = minf(l+2,l+m)
	for k in range(mink,maxk+1):
		hlp += pow(-1,l+m-k)*bino(l+2,k)*bino(l-2,k-m-2)*pow(cos(the/2),2*k-m-2)*pow(sin(the/2),2*l-2*k+m+2)

	return sqrt(fact(l+m)*fact(l-m)*(2*l+1)/(4*M_PI*fact(l-2)*fact(l+2)))*hlp*cexp(I*m*phi)


# Create the array of spherical harmonics
cdef double complex y[6][11]

cdef createy(float the, float phi):
	cdef int l, m
	
	for l in range(2,6):
		for m in range(2*l+1):
			mhlp = cmod(m,l+1)-l*(m/(l+1))
			mhlp2 = cmod(mhlp+11,11)

			y[l][mhlp2] = Y(l,mhlp,the,phi)


# Define the class as Python-code
class Fn:
	# Pass the surrogate model, the total mass (in M_sun), distance of the source (in Mpc), the sky position (in degree), the velocity magnitude (in km/s), the direction of the velocity (in degree) and the time steps (in s)
	def __init__(self, surrogate_, M_, q_, dis_, the_, phi_, v_mag_, v_the_, v_phi_, times_):
        
		self.surrogate_ = surrogate_
		self.M_ = M_*M_sun
		self.q_ = q_
		self.dis_ = dis_*Mpc
		self.the_ = pi*float(the_)/180
		self.phi_ = pi*float(phi_)/180
		self.v_mag_ = 1000*float(v_mag_)/c
		self.v_the_ = pi*float(v_the_)/180
		self.v_phi_ = pi*float(v_phi_)/180
		self.times_ = times_


	# Define the waveform
	def WaVe(self):

		sur = self.surrogate_
		M = self.M_
		q = self.q_
		dis = self.dis_
		the = self.the_
		phi = self.phi_
		v_mag = self.v_mag_
		v_the = self.v_the_
		v_phi = self.v_phi_
		times = self.times_

		times = np.array(times)
        
		# Set the spins of the BHs        
		chiA = [0,0,0]
		chiB = [0,0,0]

		# Compute velocity vector
		v = np.array([sin(v_the)*cos(v_phi),sin(v_the)*sin(v_phi),cos(v_the)])
        
		# Compute the Doppler shift and the amplitude of the wave
		Dopp = (1+v_mag*(sin(the)*cos(phi)*v[0]+sin(the)*sin(phi)*v[1]+cos(the)*v[2]))/sqrt(1-pow(v_mag,2))
		amp = G*M/(dis*pow(c,2))
        
		# Transform the times for evaluation in natural units and set f_low to zero so that it can in principle compute all times
		times = pow(c,3)*times/(G*M*Dopp)
		f_low = 0
        
		# 'Get' the wave form
		t, h, dyn = sur(q, chiA, chiB, times=times, f_low=f_low)
        
		# Compute and write out the spherical harmonics
		createy(the,phi)
		
		# Calling the class for Coeff
		Coeff_ = Coeff.Fn2(v_mag, v)
		c0, cp, cm = Coeff_.Coefficient()
        
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
						hlp0 += pow(-0.25,l0)*conj(h[(l0,abs(mi[1]))])*c0[l0][mi[0]][mi[1]]
					else:
						hlp0 += pow(0.25,l0)*h[(l0,mi[1])]*c0[l0][mi[0]][mi[1]]

			hlpp = np.zeros(len(t), dtype=complex)
			for lp in range(max(2,abs(mi[1]+1)),6):
				try:
					h[(lp,abs(mi[1]+1))]
				except KeyError:
					pass
				else:
					if mi[1]+1 < 0:
						hlpp += pow(-0.25,lp)*conj(h[(lp,abs(mi[1]+1))])*cp[lp][mi[0]][mi[1]]
					else:
						hlpp += pow(0.25,lp)*h[(lp,mi[1]+1)]*cp[lp][mi[0]][mi[1]]

			hlpm = np.zeros(len(t), dtype=complex)
			for lm in range(max(2,abs(mi[1]-1)),6):
				try:
					h[(lm,abs(mi[1]-1))]
				except KeyError:
					pass
				else:
					if mi[1]-1 < 0:
						hlpm += pow(-0.25,lm)*conj(h[(lm,abs(mi[1]-1))])*cp[lp][mi[0]][mi[1]]
					else:
						hlpm += pow(0.25,lm)*h[(lm,mi[1]-1)]*cm[lm][mi[0]][mi[1]]

			# Sum up the contributions to get the exited mode
			if mi[1] < 0:
				hb[mi] = pow(-1,mi[0])*conj(h[(mi[0],abs(mi[1]))]) + pow(-1,mi[1])*(1j/pow(4,mi[0]+1))*(hlp0+hlpp+hlpm)
			else:
				hb[mi] = h[mi] + pow(-1,mi[1])*(1j/pow(4,mi[0]+1))*(hlp0+hlpp+hlpm)

		# Define arrays for the polarizations
		hp = np.zeros(len(t), dtype=float)
		hc = np.zeros(len(t), dtype=float)
        
		# Compute polarizations in 'SI' units
		hlp = np.zeros(len(t), dtype=complex)
		for mi in hb:
			mhlp = cmod(mi[1]+11,11)

			hlp += hb[mi]*y[mi[0]][mhlp]

		hp, hc = amp*hlp.real, amp*hlp.imag

		# Write out the waveforms
		return (hp,hc)

