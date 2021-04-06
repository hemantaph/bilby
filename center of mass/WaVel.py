### Computes the waveform for moving sources based on GWSurrogate  ###

import argparse
import sys
import subprocess

import gwsurrogate
import numpy as np

from math import cos
from math import sin
from math import sqrt
from math import factorial
from math import pi
from math import e
from numpy import conj


# Get the arguments of the wave: total mass (in M_sun), mass ratio, luminosity distance (in Mpc), the sky location of the observer (in rad), the magnitude of the velocity (in km/s), the direction of the velocity (in rad), the time array (in s)
parser = argparse.ArgumentParser()
parser.add_argument('M')
parser.add_argument('q')
parser.add_argument('dis')
parser.add_argument('the')
parser.add_argument('phi')
parser.add_argument('v_mag')
parser.add_argument('v_the')
parser.add_argument('v_phi')
parser.add_argument('times')


# Define constants: speed of light (in m/s), gravitational constant (in m^3/(kg*s^2)), solar mass (in kg) and megaparsec (in m)
c = 2.99e8
G = 6.67e-11
M_sun = 1.99e30
Mpc = 3.09e22


# Set the total mass (in kg), mass ratio, luminosity distance in m and spins
M = float(parser.parse_args().M)*M_sun
q = float(parser.parse_args().q)
dis = float(parser.parse_args().dis)*Mpc
chiA = [0,0,0]
chiB = [0,0,0]

# Set sky location of observer (in deg)
the = pi*float(parser.parse_args().the)/180
phi = pi*float(parser.parse_args().phi)/180

# Set the magnitude (in c) and the direction of the velocity
v_mag = float(parser.parse_args().v_mag)/c
v_the = pi*float(parser.parse_args().v_the)/180
v_phi = pi*float(parser.parse_args().v_phi)/180
v = np.array([sin(v_the)*cos(v_phi),sin(v_the)*sin(v_phi),cos(v_the)])

# Set the times for evaluation (in s)
times = np.array(parser.parse_args().times.split()).astype(np.float)


# Compute the Doppler shift and the amplitude of the wave
Dopp = (1+v_mag*(sin(the)*cos(phi)*v[0]+sin(the)*sin(phi)*v[1]+cos(the)*v[2]))/sqrt(1-v_mag**2)
amp = G*M/(dis*c**2)

# Transform the times for evaluation in natural units and set f_low to zero so that it can in principle compute all times
for i in range(len(times)):
	times[i] = c**3*times[i]/(G*M*Dopp)

f_low = 0


# Define binomial function and spherical harmonics
def binomial(n,k):
	if n < 0 or k < 0 or n-k < 0:
		return 0
	else:
		return factorial(n)/(factorial(k)*factorial(n-k))

s = -2
def Y(l,m,THE,PHI):
	hlp = 0
	for k in range(max(0,m-s),l-s+1):
		if binomial(l-s,k) != 0 and binomial(l+s,k+s-m) != 0:
			hlp = hlp + (-1)**(l-k-s+m)*binomial(l-s,k)*binomial(l+s,k+s-m)*cos(THE/2)**(2*k+s-m)*sin(THE/2)**(2*l-2*k-s+m)*e**(1j*m*PHI)

	return sqrt(factorial(l+m)*factorial(l-m)*(2*l+1)/(4*pi*factorial(l+s)*factorial(l-s)))*hlp


# Load the surrogate model
sur = gwsurrogate.LoadSurrogate('NRHybSur3dq8')

# 'Get' the wave form
t, h, dyn = sur(q, chiA, chiB, times=times, f_low=f_low)


# Compute and write out the spherical harmonics
y = np.zeros((6,12), dtype=complex)

for l2 in range(2,6):
	for m in range(-l2,l2+1):
		y[l2][m] = Y(l2,m,the,phi)


# Convert 'v' into a string, run 'coeff.py' and import the velocity coefficients
v_str = ' '.join(map(str,v))
subprocess.run([sys.executable, "coeff.py", str(v_mag), v_str])
coeff = np.load('coeff.npy')
c0, cp, cm = coeff[0], coeff[1], coeff[2]
			

# Compute the beamed modes
# Create array for beamed modes
hb = {}
for mi in h:
	hb[mi] = [0]*len(t)
	hb[(mi[0],-mi[1])] = [0]*len(t)

# Compute the beamed modes over time
for mi in hb:
	for i in range(len(t)):
		# Compute the contributions of each mode
		hlp0 = 0
		for l0 in range(max(2,abs(mi[1])),6):
			try:
				h[(l0,abs(mi[1]))]
			except KeyError:
				pass
			else:
				if mi[1] < 0:
					hlp0 = hlp0 + (1/4**l0)*(-1)**l0*conj(h[(l0,abs(mi[1]))][i])*c0[l0][mi[0]][mi[1]]
				else:
					hlp0 = hlp0 + (1/4**l0)*h[(l0,mi[1])][i]*c0[l0][mi[0]][mi[1]]
				
		hlpp = 0
		for lp in range(max(2,abs(mi[1]+1)),6):
			try:
				h[(lp,abs(mi[1]+1))]
			except KeyError:
				pass
			else:
				if mi[1]+1 < 0:
					hlpp = hlpp + (1/4**lp)*(-1)**lp*conj(h[(lp,abs(mi[1]+1))][i])*cp[lp][mi[0]][mi[1]]
				else:
					hlpp = hlpp + (1/4**lp)*h[(lp,mi[1]+1)][i]*cp[lp][mi[0]][mi[1]]

		hlpm = 0
		for lm in range(max(2,abs(mi[1]-1)),6):
			try:
				h[(lm,abs(mi[1]-1))]
			except KeyError:
				pass
			else:
				if mi[1]-1 < 0:
					hlpm = hlpm + (1/4**lm)*(-1)**lm*conj(h[(lm,abs(mi[1]-1))][i])*cp[lp][mi[0]][mi[1]]
				else:
					hlpm = hlpm + (1/4**lm)*h[(lm,mi[1]-1)][i]*cm[lm][mi[0]][mi[1]]
		
		# Sum up the contributions to get the exited mode
		if mi[1] < 0:
			hb[mi][i] = (-1)**mi[0]*conj(h[(mi[0],abs(mi[1]))][i]) + (-1)**abs(mi[1])*(1j/4**(mi[0]+1))*(hlp0+hlpp+hlpm)
		else:
			hb[mi][i] = h[mi][i] + (-1)**mi[1]*(1j/4**(mi[0]+1))*(hlp0+hlpp+hlpm)


# Define arrays for the waveforms
hp = [0]*len(t)
hc = [0]*len(t)

# Compute polarizations over time
for i in range(len(t)):
	
	hlp = 0
	for mi in hb:
		hlp = hlp + hb[mi][i]*y[mi[0]][mi[1]]

	hp[i], hc[i] = hlp.real, hlp.imag

	# Convert time and amplitude in 'SI' units
	hp[i], hc[i] = amp*hp[i], amp*hc[i]


# Write out the waveforms	
np.save('WaVel',[hp,hc])

