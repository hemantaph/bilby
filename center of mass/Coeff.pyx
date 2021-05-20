### Computes coefficients for the modes excited by motion  ###

import numpy as np

from cython import cmod

from libc.math cimport sqrt, pow, abs
from libc.math cimport M_PI

cdef extern from "complex.h":
	float complex I


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


# Define the class as Python-code
class Fn2:
	def __init__(self, v_mag_, v_):
        
		self.v_mag_ = v_mag_
		self.v_ = v_

	# Compute coefficients und save them in an array
	def Coefficient(self):
		v_mag = self.v_mag_
		v = self.v_

		# Compute the coefficients
		createc(v_mag,v[0],v[1],v[2])

		return (coeff[0],coeff[1],coeff[2])


