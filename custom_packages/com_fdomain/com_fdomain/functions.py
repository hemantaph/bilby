from bilby.gw.source import lal_binary_black_hole
from bilby.gw.waveform_generator import WaveformGenerator
import gwsurrogate
import gwpy
import numpy as np

from com_fdomain.cython_files import WaVel

surrogate = gwsurrogate.LoadSurrogate('NRHybSur3dq8')

#hard coding of time series
#becasue the input is frequency domain but we need time series to find h+x 
#Also the timeseries is calculated outside the function inorder to avoid repetition in each iteration
#dt = 1/sampling_frequency
#time_array_ = np.arange(start_time,start_time+duration,dt)
#start_time = 0.0-duration+(1/4*duration)
#duration=4s
dt = 1/1024
time_array_ = np.arange(-3.98,0.02,dt)

#####################################################################
# Define the frequency-domain model
#geocent_time=1126259600.
#trigger-time=1126259598.02
#start_time_gps = geocent_time-4.0+0.02
start_time_gps = 126259600.0-4.0+0.02
def moving_bbh(freq_array, total_mass, mass_ratio, luminosity_distance, theta_jn, psi, speed, v_the, v_phi,**kwargs):
    #len(times)=duration*sampling_frequency
    times = time_array_
    
    arg_ = {'surrogate_':surrogate,'M_':total_mass, 'q_':mass_ratio, 'dis_':luminosity_distance, 'the_':theta_jn, 'phi_':psi, 'v_mag_':speed, 'v_the_':v_the, 'v_phi_':v_phi, 'times_':times, **kwargs}
    fplus = WaVel.Fn(**arg_)
    h_plus = gwpy.timeseries.TimeSeries(fplus.WaVe()[0], t0=start_time_gps,\
                                        dt=dt, name='Strain', channel=None)
    h_cross = gwpy.timeseries.TimeSeries(fplus.WaVe()[1], t0=start_time_gps,\
                                        dt=dt, name='Strain', channel=None)
    
    fplus = h_plus.fft()
    fcross = h_cross.fft()

    return {'plus': fplus.value, 'cross': fcross.value}
  
#####################################################################
#waveform-generator with frequency_domain_source_model
#the argugent values are taken from inside bilby pipe
def my_waveform_generator(duration, sampling_frequency, frequency_domain_source_model, **kwargs):
    
    return WaveformGenerator(duration=duration, sampling_frequency=sampling_frequency, frequency_domain_source_model=frequency_domain_source_model)

#####################################################################