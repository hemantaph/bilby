from bilby.gw.source import lal_binary_black_hole
from bilby.gw.waveform_generator import WaveformGenerator
import numpy as np
#from bilby.gw.likelihood import GravitationalWaveTransient

def my_function(time, amplitude, damping_time, frequency, phase, t0):
    plus = np.zeros(len(time))
    tidx = time >= t0
    plus[tidx] = amplitude * np.exp(-(time[tidx] - t0) / damping_time) *\
        np.sin(2 * np.pi * frequency * (time[tidx] - t0) + phase)
    cross = np.zeros(len(time))
    return {'plus': plus, 'cross': cross}


def my_waveform_generator(duration, sampling_frequency, frequency_domain_source_model, **kwargs):
    return WaveformGenerator(duration=duration, sampling_frequency=sampling_frequency, time_domain_source_model=time_domain_source_model, start_time=1126259600. - 0.5)
