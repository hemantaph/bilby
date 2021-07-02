from bilby.gw.waveform_generator import WaveformGenerator
import numpy as np
import eccn.hphc as hphc

C = 299792458.
G = 6.67408*1e-11
Mo = 1.989*1e30
Mpc = 3.086*1e22
#####################################################################
# Define the frequency-domain waveform model
def eccentric_waveform(frequency_array_, chirp_mass, mass_ratio, eccentricity, luminosity_distance, theta_jn, psi, phase, geocent_time, ra, dec, **kwargs):
    
    mass_1 = (chirp_mass*(1+mass_ratio)**(1/5))/mass_ratio**(3/5)
    mass_2 = chirp_mass*mass_ratio**(2/5)*(1+mass_ratio)**(1/5)
    luminosity_distance = luminosity_distance*Mpc
    total_mass = (mass_1+mass_2)*Mo
    symmetric_mass_ratio = (mass_1*mass_2)/((mass_1+mass_2)**2)
    lso_f = (C**3)/( G*(mass_1+mass_2)*Mo*np.pi*6**(3/2) )
    mass_diff = (mass_1-mass_2)*Mo
    f_min = 20.
    f_max = lso_f
    
    foo = np.array(frequency_array_, dtype='float')

    #h_plus,h_cross = htilde(theta_jn, psi, luminosity_distance , foo, f_min, eccentricity, phase, luminosity_distance, total_mass, symmetric_mass_ratio, lso_f, mass_diff)
    
    arg_plus = {'iota_':theta_jn, 'beta_':psi, 'D_':luminosity_distance , \
                'farray_':foo, 'f0_':20.0, 'et0_':eccentricity, 'phic_':phase, \
                'tc_':geocent_time, 'M_':total_mass, 'eta_':symmetric_mass_ratio, \
                'ff_':lso_f, 'delta_':mass_diff}

    fplus = hphc.Fn(**arg_plus)

    h_plus = fplus.htilde()[0]
    h_cross = fplus.htilde()[1]

    return {'plus': h_plus, 'cross': h_cross}
    #return(hphc.htilde(theta_jn, psi, luminosity_distance , foo, f_min, eccentricity, phase, geocent_time, total_mass, symmetric_mass_ratio, lso_f, mass_diff, f_max))
#####################################################################

def my_waveform_generator(duration, sampling_frequency, frequency_domain_source_model, **kwargs):
    
    return WaveformGenerator(duration=duration, sampling_frequency=sampling_frequency, frequency_domain_source_model=frequency_domain_source_model, **kwargs)

#####################################################################