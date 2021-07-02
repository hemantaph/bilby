from bilby.gw.source import lal_binary_black_hole
from bilby.gw.waveform_generator import WaveformGenerator
import gwsurrogate
#import gwpy
#import numpy

import WaVel

sampling_frequency = 1024
duration = 2.0
geocent_time = 1126259462.4
start_time = geocent_time - duration + 0.02

surrogate = gwsurrogate.LoadSurrogate('NRHybSur3dq8')

#####################################################################
# Define the time-domain model
def moving_bbh(times, total_mass, mass_ratio, luminosity_distance, theta, phi, speed, v_the, v_phi,**kwargs):

    arg_ = {'surrogate_':surrogate,'M_':total_mass, 'q_':mass_ratio, 'dis_':luminosity_distance, 'the_':theta, 'phi_':phi, 'v_mag_':speed, 'v_the_':v_the, 'v_phi_':v_phi, 'times_':times, **kwargs}
    
    fplus = WaVel.Fn(**arg_)
    h_plus = fplus.WaVe()[0]
    h_cross = fplus.WaVe()[1]
    
    return {'plus': h_plus, 'cross': h_cross}
#####################################################################

def my_waveform_generator(duration, sampling_frequency, **kwargs):
    
    return WaveformGenerator(duration=duration, sampling_frequency=sampling_frequency, time_domain_source_model=moving_bbh, start_time=start_time)

#####################################################################











