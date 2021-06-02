from bilby.gw.source import lal_binary_black_hole
from bilby.gw.waveform_generator import WaveformGenerator
#from bilby.gw.likelihood import GravitationalWaveTransient


def my_function(frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
                phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, **kwargs):
    return lal_binary_black_hole(frequency_array, mass_1, mass_2, luminosity_distance,
                                 a_1, tilt_1, phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, **kwargs)

def my_waveform_generator(duration, sampling_frequency, **kwargs):
    
    return WaveformGenerator(duration=duration, sampling_frequency=sampling_frequency, frequency_domain_source_model=my_function)

'''
def my_likelihood(duration, sampling_frequency, frequency_domain_source_model, **kwargs):
    
    return GravitationalWaveTransient(duration=duration, sampling_frequency=sampling_frequency, frequency_domain_source_model=frequency_domain_source_model)
'''