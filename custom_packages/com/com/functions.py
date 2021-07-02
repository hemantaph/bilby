from bilby.gw.source import lal_binary_black_hole
from bilby.gw.waveform_generator import WaveformGenerator
import gwsurrogate
#import gwpy
#import numpy

from com.c_files import WaVel


geocent_time_ = 1126259462.4

surrogate = gwsurrogate.LoadSurrogate('NRHybSur3dq8')

#####################################################################
# Define the time-domain model
def moving_bbh(times, mass, ratio, distance, theta, phi, speed, v_the, v_phi,**kwargs):
    
    times = times - geocent_time_
    arg_ = {'surrogate_':surrogate,'M_':mass, 'q_':ratio, 'dis_':distance, 'the_':theta, 'phi_':phi, 'v_mag_':speed, 'v_the_':v_the, 'v_phi_':v_phi, 'times_':times, **kwargs}
    fplus = WaVel.Fn(**arg_)
    h_plus = fplus.WaVe()[0]
    h_cross = fplus.WaVe()[1]
   
    return {'plus': h_plus, 'cross': h_cross}
    
#####################################################################

def my_waveform_generator(duration, sampling_frequency, **kwargs):
    
    return WaveformGenerator(duration=duration, sampling_frequency=sampling_frequency, time_domain_source_model=moving_bbh)

#####################################################################