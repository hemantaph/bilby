#create_your_own_source_model

from bilby.gw.source import lal_binary_black_hole
from bilby.gw.waveform_generator import WaveformGenerator
import numpy as np
#from bilby.gw.likelihood import GravitationalWaveTransient

'''
def my_function(frequency_array, A, f0, tau, phi0, geocent_time, ra, dec, psi):
    
    arg = -(np.pi * tau * (f - f0))**2 + 1j * phi0
    plus = np.sqrt(np.pi) * A * tau * np.exp(arg) / 2.
    cross = plus * np.exp(1j * np.pi / 2)
    return {'plus': plus, 'cross': cross}


def my_function(frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
                phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, **kwargs):
    n = len(frequency_array)
    plus = np.zeros(n)*1j
    cross = np.zeros(n)*1j
    return {'plus': plus, 'cross': cross}


def my_function(frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
                phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, **kwargs):
    return lal_binary_black_hole(frequency_array, mass_1, mass_2, luminosity_distance,
                                 a_1, tilt_1, phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, **kwargs)


def my_function(frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
        phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, A, f0, tau, phi0):
    
    arg = -(np.pi * tau * (f - f0))**2 + 1j * phi0
    plus = np.sqrt(np.pi) * A * tau * np.exp(arg) / 2.
    cross = plus * np.exp(1j * np.pi / 2)
    return {'plus': plus, 'cross': cross}
'''
def my_function(frequency_array, A, f0, tau, phi0, geocent_time, ra, dec, psi, **kwargs):
    
    f = frequency_array
    A = A*1e-23
    arg = -(np.pi * tau * (f - f0))**2 + 1j * phi0
    plus = np.sqrt(np.pi) * A * tau * np.exp(arg) / 2.
    cross = plus * np.exp(1j * np.pi / 2)
    return {'plus': plus, 'cross': cross}

def my_waveform_generator(duration, sampling_frequency, frequency_domain_source_model, **kwargs):
    
    return WaveformGenerator(duration=duration, sampling_frequency=sampling_frequency, frequency_domain_source_model=frequency_domain_source_model)


