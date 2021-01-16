#!/usr/bin/env python
"""
A script to show how to create your own time domain source model.
A simple damped Gaussian signal is defined in the time domain, injected into
noise in two interferometers (LIGO Livingston and Hanford at design
sensitivity), and then recovered.
"""
from __future__ import division, print_function

import numpy as np
import bilby
import lal
from scipy.signal import tukey
from pycbc.waveform import ringdown as rd



event = 'GW150914'
time_of_event = bilby.gw.utils.get_event_time(event)
print("Time of Event= "+str(time_of_event))

#time_of_event = 1126259462

# This sets up logging output to understand what bilby is doing
#bilby.core.utils.setup_logger(outdir=outdir, label=label)

# Here we import the detector data. This step downloads data from the
# LIGO/Virgo open data archives. The data is saved to an `Interferometer`
# object (here called `H1` and `L1`). A Power Spectral Density (PSD) estimate
# is also generated and saved to the same object. To check the data and PSD
# makes sense, for each detector a plot is created in the `outdir` called
# H1_frequency_domain_data.png and LI_frequency_domain_data.png. The two
# objects are then placed into a list.
# For GW170817, 170608 and 170814 add the following line to select clean data.
# kwargs = {"tag": 'CLN'}




# define the time-domain model

def Ydirect(i,p,l,m):
    Y_lm = lal.SpinWeightedSphericalHarmonic(i, p, -2, l, m)
    return Y_lm

def ringdown(
        time, final_mass, final_spin, iota, varphi, amplitude,phase,amplitude_1,phase_1,amplitude_2,phase_2,amplitude_3,phase_3,amplitude_4,phase_4,amplitude_5,phase_5,amplitude_6,phase_6,amplitude_7,phase_7,geocent_time):
    """
    This example only creates a linearly polarised signal with only plus
    polarisation.i
    """

    f,tau = rd.get_lm_f0tau_allmodes(final_mass,final_spin,['228'])

    frequency_220 = f['220']
    tau_220 = tau['220']
    frequency_221 = f['221']
    tau_221 = tau['221']
    frequency_222 = f['222']
    tau_222 = tau['222']
    frequency_223 = f['223']
    tau_223 = tau['223']
    frequency_224 = f['224']
    tau_224 = tau['224']
    frequency_225 = f['225']
    tau_225 = tau['225']
    frequency_226 = f['226']
    tau_226 = tau['226']
    frequency_227 = f['227']
    tau_227 = tau['227']




    h22 = (amplitude) * np.exp(-(np.abs(time-geocent_time)) / tau_220) *\
    np.exp(1j *(2 * np.pi * frequency_220 * time + phase)) + (amplitude_1*amplitude) * np.exp(-(np.abs(time-geocent_time)) / tau_221) *\
    np.exp(1j *(2 * np.pi * frequency_221 * time + phase_1)) + (amplitude_2) * np.exp(-(np.abs(time-geocent_time)) / tau_222) *\
    np.exp(1j *(2 * np.pi * frequency_222 * time + phase_2)) + (amplitude_3) * np.exp(-(np.abs(time-geocent_time)) / tau_223) *\
    np.exp(1j *(2 * np.pi * frequency_223 * time + phase_3)) + (amplitude_4) * np.exp(-(np.abs(time-geocent_time)) / tau_224) *\
    np.exp(1j *(2 * np.pi * frequency_224 * time + phase_4)) + (amplitude_5) * np.exp(-(np.abs(time-geocent_time)) / tau_225) *\
    np.exp(1j *(2 * np.pi * frequency_225 * time + phase_5)) + (amplitude_6) * np.exp(-(np.abs(time-geocent_time)) / tau_226) *\
    np.exp(1j *(2 * np.pi * frequency_226 * time + phase_6)) + (amplitude_7) * np.exp(-(np.abs(time-geocent_time)) / tau_227) *\
    np.exp(1j *(2 * np.pi * frequency_227 * time + phase_7))

    h2m2 = np.conj(h22)


    h = Ydirect(iota,varphi,2,2) * h22 + Ydirect(iota,varphi,2,-2) * h2m2


    #cross =  np.zeros(len(time)) 
    plus = np.real(h)
    cross = np.imag(h)


    template_duration = 0.125
    start_index = np.argmin(np.abs(time-geocent_time))
    end_index = np.argmin(np.abs(time-geocent_time-template_duration))

    for i in np.arange(0,start_index+1):
        plus[i]=0
        cross[i]=0


    #for i in np.arange(end_index,len(time)):
    #    plus[i]=0
    #    cross[i]=0

    window_alpha = 0.00/template_duration
    window_start_point = 0.

    window = tukey(end_index-start_index+1,alpha=2*window_alpha)
    mid_index = np.argmin(np.abs(time-geocent_time-window_start_point))

    #print('Len of template window= '+str(len(window)))

    for i in np.arange(start_index,end_index+1):

             if i < mid_index:

                      plus *= 1
             else:

                      plus[i] *= window[i - start_index]


    #window = tukey(len(time),alpha=0.)

    
    return {'plus': plus, 'cross': cross}
 



def time_domain_damped_sinusoid_no_hair(
        time, amplitude, damping_time, frequency, phase,amplitude_2,frequency_2, damping_time_2, phase_2,geocent_time):
    """
    This example only creates a linearly polarised signal with only plus
    polarisation.i
    """
   
    if amplitude_2 > amplitude:

 	plus = np.zeros(len(time))
        

    else:
    

    	plus = (amplitude) * np.exp(-(np.abs(time+0.012)) / damping_time) *\
    	np.sin(2 * np.pi * frequency * time + phase) + (amplitude_2) * np.exp(-(np.abs(time+0.012)) /(damping_time_2)) *\
    	np.sin(2 * np.pi * frequency_2 * time + phase_2)
    
    
    #cross =  np.zeros(len(time)) 
        
    window = tukey(len(time),alpha=0.)
    
    plus *= window
    cross = np.zeros(len(time))

    #cross = plus
    return {'plus': plus, 'cross': cross}


def time_domain_damped_sinusoid_single_mode(
        time, amplitude, damping_time, frequency, phase):
    """
    This example only creates a linearly polarised signal with only plus
i    polarisation.
    """

   

    plus = amplitude * np.exp(-(time) / damping_time) *\
    np.sin(2 * np.pi * frequency * time + phase) 


    cross = np.zeros(len(time))
    return {'plus': plus, 'cross': cross}

pre_peak_time = 0
duration = 4
sampling_frequency = 4096
outdir = '/home/juancalderonbustillo/public_html/NoHair/PreProd_5.4'
label = 'Kerr_2mode_ISISkyLoc_Prod_1024_4s_PSD_CheckGreg_500s'
NR_path = '/home/gayathri.v/git_gayathri/lvcnr-lfs/SXS/SXS_BBH_0305_Res5.h5'

distance = 0.7*400/3
#distance=0.091*400
#distance = 30/0.091
# define parameters to inject.

injection_parameters = dict(
mtotal=72., luminosity_distance=distance, iota=0, psi=0.,
    phase=1.3, geocent_time=0., ra=1.375, dec=-1.2108)

waveform_arguments = dict(waveform_approximant='NRhdf5',
                          reference_frequency=0., minimum_frequency=20,t_start=0.00,alpha=0.0,nr_file = NR_path,selected_modes = [[2,2]])

waveform_injection_generator = bilby.gw.WaveformGenerator(
    duration=duration, sampling_frequency=sampling_frequency,
    time_domain_source_model=bilby.gw.source.pycbc_nr_binary_black_hole_time_domain,
    waveform_arguments=waveform_arguments,start_time=time_of_event-pre_peak_time)


#template_fake_parameters = dict(amplitude=1e-21, damping_time=0.01, frequency=200,
                            #phase=0, geocent_time=injection_parameters['geocent_time'],ra=1.375, dec=-1.2108,
			#	psi=0)


template_fake_parameters = dict(final_mass=69, final_spin=0.67,iota=np.pi,varphi=0.,amplitude = 1e-21,
                            phase=0,amplitude_1 = 0,
                            phase_1=0,amplitude_2 = 0,
                            phase_2=0,amplitude_3 = 0,
                            phase_3=0,amplitude_4 = 0,
                            phase_4=0,amplitude_5 = 0,
                            phase_5=0,amplitude_6 = 0,
                            phase_6=0,amplitude_7 = 0,
                            phase_7=0,geocent_time=1126259462.4176,#ra=1.375, dec=-1.2108,psi=0,
			    ra=1.95, dec=-1.27,psi=0.82)


#ra=1.375, dec=-1.2108,psi=0


waveform = bilby.gw.waveform_generator.WaveformGenerator(
    duration=duration, sampling_frequency=sampling_frequency,
    time_domain_source_model=ringdown,start_time=time_of_event-pre_peak_time)


ringdown_waveform = bilby.gw.waveform_generator.WaveformGenerator(
    duration=duration, sampling_frequency=sampling_frequency,
    time_domain_source_model=ringdown,start_time=time_of_event-0.9
)


# GET DATA FROM INTERFEROMETER
interferometer_names = ['H1','L1']  # can also include 'V1'
#duration = 4    # length of data segment containing the signal
roll_off = 0.2  # smoothness of transition from no signal
# to max signal in a Tukey Window.
psd_offset = -1024  # PSD is estimated using data from
# 'center_time + psd_offset' to 'center_time + psd_offset + psd_duration'
psd_duration = 500
#filter_freq = None  # low pass filter frequency to cut signal content above
# Nyquist frequency. The condition is 2 * filter_freq >= sampling_frequency

# Keyword args are passed to 'gwpy.timeseries.TimeSeries.fetch_open_data()'
# sample_rate: most data are stored by LOSC at 4096 Hz,
# there may be event-related data releases with a 16384 Hz rate.
kwargs = {"sample_rate":4096}
# For some events a "tag" is required to download the data. (CLN = clean data)
#kwargs = {"tag": 'CLN'}
# For some events can specify channels: source data stream for LOSC data.
#kwargs = {"channel": {'H1': 'H1:DCS-CALIB_STRAIN_C02'}}
#                      'L1': 'L1:DCS-CALIB_STRAIN_C02',
#                      #'V1': 'V1:FAKE_h_16384Hz_4R'
#		}}

#option event removed

#ifos = bilby.gw.detector.get_interferometer_with_open_data('H1',time_of_event,duration=duration, start_time=time_of_event, roll_off=roll_off, outdir=outdir, label=label)


ifos = bilby.gw.detector.get_event_data(event,
    #start_time = time_of_event-0.01,
    start_time = 1126259462.4-2.0,
    interferometer_names=interferometer_names,
    duration=duration,
    roll_off = roll_off,
    psd_duration=psd_duration,plot=False, psd_offset = psd_offset #,minimum_frequency = 60
    #filter_freq = None,**kwargs
    )

#ifos.power_spectral_density = bilby.gw.detector.PowerSpectralDensity(psd_file='./outdir/H1_PSD_1126258439.5_100.0.txt')

ifos.plot_data(outdir='/home/juancalderonbustillo/public_html/NoHair')

# call the waveform_generator to create our waveform model.

# inject the signal into three interferometers
#ifos = bilby.gw.detector.InterferometerList(['H1'])

#for ifo in ifos:
#	ifo.set_strain_data_from_zero_noise(
#    	sampling_frequency=sampling_frequency, duration=duration,
# 	   start_time=injection_parameters['geocent_time']-pre_peak_time)

#    ifos.set_strain_data_from_power_spectral_densities(
#    sampling_frequency=sampling_frequency, duration=duration,
#    start_time=injection_parameters['geocent_time'] - pre_peak_time)

#ifos.inject_signal(waveform_generator=waveform_injection_generator,
#                   	parameters=injection_parameters)
#


#print("Plotting Injection and one template before injecting in noise")
import numpy as np
import matplotlib.pyplot as plt
plt.close()
#wf = waveform_injection_generator.time_domain_strain(injection_parameters)
ds = ringdown_waveform.time_domain_strain(template_fake_parameters)
#plt.plot(waveform_injection_generator.time_array, wf['plus'],color="red")
plt.plot(ringdown_waveform.time_array,ds['plus'],color='red')
plt.xlim(time_of_event-0.01,time_of_event+0.125)
plt.axvline(x=time_of_event,color='red')
plt.savefig('/home/juancalderonbustillo/public_html/NoHair/GW150914_ringdown_test_waveform_test.png')
plt.close()

#for ifo in ifos:
#	ifo.plot_time_domain_data(outdir=outdir)
#	print ifo.strain_data.keys()
#  create the priors
prior = template_fake_parameters.copy()
#prior['amplitude'] = bilby.core.prior.Uniform(1e-22, 1e-18, r'$h_0$')
#prior['damping_time'] = bilby.core.prior.Uniform(0.001, 0.1, r'damping time', unit='$s$')
#prior['damping_time_2'] = bilby.core.prior.Uniform(0.001, 0.1, r'damping time 2', unit='$s$')
#prior['frequency'] = bilby.core.prior.Uniform(100,500, r'frequency', unit='Hz')
#prior['frequency_2'] = bilby.core.prior.Uniform(100,500, r'frequency_2', unit='Hz')
#prior['amplitude_2'] = bilby.core.prior.Uniform(1e-22, 1e-19, r'$h2_0$')
#prior['phase_2'] = bilby.core.prior.Uniform(0, 2*np.pi, r'$\phi_2$')
#prior['phase'] = bilby.core.prior.Uniform(0, 2*np.pi, r'$\phi$')
#prior['psi'] = bilby.core.prior.Uniform(-np.pi / 2, np.pi / 2, r'$\psi$')
#prior['ra']=bilby.core.prior.Uniform(0, 2*np.pi, r'$ra$')
#prior['dec']=bilby.core.prior.Uniform(-np.pi / 2, np.pi / 2, r'$dec$')
#prior['geocent_time'] = bilby.core.prior.Uniform(injection_parameters['geocent_time']-0.01, injection_parameters['geocent_time']+0.01, r'$t_c$')


#Mass and spin priors for ringdown analysis

prior['final_mass'] = bilby.core.prior.Uniform(20, 100, r'$M_{f}$')
prior['final_spin'] = bilby.core.prior.Uniform(0., 0.98, r'$a_{f}$', unit='$s$')
#prior['geocent_time'] = bilby.core.prior.Uniform(1126259462.42 ,1126259462.44, r'$t_{c}$', unit='$s$')
prior['geocent_time'] = bilby.core.prior.Uniform(time_of_event , time_of_event+0.4, r'$t_{c}$', unit='$s$')
#prior['geocent_time'] = time_of_event
#Damped Sinusiod priors
#prior['iota'] = bilby.core.prior.Sine(name='iota')
prior['varphi'] = bilby.core.prior.Uniform(0, 2*np.pi, r'$\varphi$')


prior['amplitude'] = bilby.core.prior.Uniform(1e-22, 1e-18, r'$h_0$')
#prior['amplitude_x'] = bilby.core.prior.Uniform(0, 1, r'$h_X$')
prior['phase'] = bilby.core.prior.Uniform(0, 2*np.pi, r'$\phi$')
prior['amplitude_1'] = bilby.core.prior.Uniform(0.1, 100, r'$h_01$')
prior['phase_1'] = bilby.core.prior.Uniform(0, 2*np.pi, r'$\phi_1$')
#prior['amplitude_2'] = bilby.core.prior.Uniform(1e-24, 1e-18, r'$h_02$')
#prior['phase_2'] = bilby.core.prior.Uniform(0, 2*np.pi, r'$\phi_2$')
#prior['amplitude_3'] = bilby.core.prior.Uniform(1e-24, 1e-18, r'$h_03$')
#prior['phase_3'] = bilby.core.prior.Uniform(0, 2*np.pi, r'$\phi_3$')
#prior['amplitude_4'] = bilby.core.prior.LogUniform(1e-24, 1e-18, r'$h_04$')
#prior['phase_4'] = bilby.core.prior.Uniform(0, 2*np.pi, r'$\phi_4$')
#prior['amplitude_5'] = bilby.core.prior.LogUniform(1e-24, 1e-18, r'$h_05$')
#prior['phase_5'] = bilby.core.prior.Uniform(0, 2*np.pi, r'$\phi_5$')
#prior['amplitude_6'] = bilby.core.prior.LogUniform(1e-24, 1e-18, r'$h_06$')
#prior['phase_6'] = bilby.core.prior.Uniform(0, 2*np.pi, r'$\phi_6$')
#prior['amplitude_7'] = bilby.core.prior.LogUniform(1e-24, 1e-18, r'$h_07$')
#prior['phase_7'] = bilby.core.prior.Uniform(0, 2*np.pi, r'$\phi_7$')

# define likelihood
likelihood = bilby.gw.likelihood.GravitationalWaveTransient(ifos, ringdown_waveform ,time_marginalization=True,priors=prior
		)

# launch sampler
result = bilby.core.sampler.run_sampler(
    likelihood, prior, sampler='cpnest', #npoints=100,
    maxmcmc=1024,outdir=outdir,nlive=1024, label=label,n_periodic_checkpoint=500
)
result.plot_corner()


