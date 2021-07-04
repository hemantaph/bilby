from bilby.gw.waveform_generator import WaveformGenerator
import numpy as np
import eccn.hphc as hphc

C = 299792458.
G = 6.67408*1e-11
Mo = 1.989*1e30
Mpc = 3.086*1e22
#####################################################################
# Define the frequency-domain waveform model
def eccentric_waveform(frequency_array_, chirp_mass, mass_ratio, initial_eccentricity, luminosity_distance, theta_jn, psi, phase, geocent_time, ra, dec, **kwargs):
    
    mass_1 = (chirp_mass*(1+mass_ratio)**(1/5))/mass_ratio**(3/5)
    mass_2 = chirp_mass*mass_ratio**(2/5)*(1+mass_ratio)**(1/5)
    luminosity_distance = luminosity_distance*Mpc
    total_mass = (mass_1+mass_2)*Mo
    symmetric_mass_ratio = (mass_1*mass_2)/((mass_1+mass_2)**2)
    #last stable orbit
    lso_f = (C**3)/( G*(mass_1+mass_2)*Mo*np.pi*6**(3/2) )
    mass_diff = (mass_1-mass_2)*Mo
    f_min = 20.
    #f_max is set according to the prior we choose, depends on total_mass. 
    f_max = 110.
    #time of coalescence is taken to be 0 for convenience
    tc = 0.
    
    foo = np.array(frequency_array_, dtype='float')

    h_plus,h_cross = hphc.htilde(foo, total_mass, symmetric_mass_ratio, mass_diff, initial_eccentricity, luminosity_distance, theta_jn, psi, phase, tc, f_min, lso_f, f_max)

    return {'plus': h_plus, 'cross': h_cross}
#####################################################################

def my_waveform_generator(duration, sampling_frequency, frequency_domain_source_model, **kwargs):
    
    return WaveformGenerator(duration=duration, sampling_frequency=sampling_frequency, frequency_domain_source_model=frequency_domain_source_model, **kwargs)

#####################################################################















