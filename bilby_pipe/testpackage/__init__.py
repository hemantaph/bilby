import numpy as np 

class test:
    def __init__(self):
        print("constructor made!!")
    
    def prt(self, number)
        print( np.power(number,2) )
        
    # frequency domain.
    def sine_gaussian(self, f, A, f0, tau, phi0, geocent_time, ra, dec, psi):
        arg = -(np.pi * tau * (f - f0))**2 + 1j * phi0
        plus = np.sqrt(np.pi) * A * tau * np.exp(arg) / 2.
        cross = plus * np.exp(1j * np.pi / 2)
        return {'plus': plus, 'cross': cross}