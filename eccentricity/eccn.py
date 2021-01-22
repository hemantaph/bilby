import numpy as np
import bilby
import hphc

C = 299792458.
G = 6.67408*1e-11
Mo = 1.989*1e30
Mpc = 3.086*1e22

duration = 2.0
sampling_frequency = 1024.0

outdir = 'outdir'
label = 'eccentric'
bilby.core.utils.setup_logger(outdir=outdir, label=label)

np.random.seed(1234)

def eccentric_waveform(frequency_array_, mass_1, mass_2, eccentricity, luminosity_distance, theta_jn, psi, phase, geocent_time, ra, dec):
    
    minimum_f = 20.0
    maximum_f = (C**3)/( G*(mass_1+mass_2)*Mo*np.pi*6**(3/2) )
    
    frequency_array = frequency_array_[(frequency_array_>=minimum_f) & (frequency_array_<=maximum_f)]
    before = frequency_array_[(frequency_array_<minimum_f)]*0.0
    after = frequency_array_[(frequency_array_>maximum_f)]*0.0
    
    N = len(frequency_array)
    h_plus = np.zeros(N)*1j
    h_cross = np.zeros(N)*1j
    
    k = 0
    for f in frequency_array:

        arg_plus = {'iota_':theta_jn, 'beta_':psi, 'D_':luminosity_distance , 'm1_':mass_1, 'm2_':mass_2, 'f_':f, 'f0_':20.0, 'et0_':eccentricity, 'phic_':phase, 'tc_':geocent_time}

        fplus = hphc.Fn(**arg_plus)

        h_plus[k] = fplus.htilde()[0]
        h_cross[k] = fplus.htilde()[1]

        k=k+1
    #return {'plus': h_plus, 'cross': h_cross}
    return {'plus': np.concatenate((before,h_plus,after)), 'cross': np.concatenate((before,h_cross,after))}




injection_parameters = dict(mass_1=36.0, mass_2=29.0, eccentricity=0.1, luminosity_distance=4000.0, theta_jn=0.4, psi=2.659, phase=1.3, geocent_time=1126259642.413, ra=1.375, dec=-1.2108)


waveform_generator = bilby.gw.waveform_generator.WaveformGenerator(
    duration=duration, sampling_frequency=sampling_frequency,
    frequency_domain_source_model=eccentric_waveform)


# Set up interferometers.
minimum_frequency = 20.0+5.0
MM1 = 35.0
MM2 = 30.0
maximum_frequency = (C**3)/( G*(MM1+MM2)*Mo*np.pi*6**(3/2) )-5.0

ifos = bilby.gw.detector.InterferometerList(['H1', 'L1'])
for ifo in ifos:
    ifo.minimum_frequency = minimum_frequency
    ifo.maximum_frequency = maximum_frequency
ifos.set_strain_data_from_power_spectral_densities(
    sampling_frequency=sampling_frequency, duration=duration,
    start_time=injection_parameters['geocent_time'] - 3)
injection = ifos.inject_signal(waveform_generator=waveform_generator,
                   parameters=injection_parameters)

# Now we set up the priors on each of the binary parameters.
prior = injection_parameters.copy()
prior['eccentricity'] = bilby.core.prior.Uniform(0.01, 0.2, 'e')

# Initialising the likelihood function.
likelihood = bilby.gw.likelihood.GravitationalWaveTransient(
    interferometers=ifos, waveform_generator=waveform_generator)

# Now we run sampler (PyMultiNest in our case).
result = bilby.core.sampler.run_sampler(
    likelihood, prior, sampler='dynesty', outdir=outdir, label=label,
    resume=False, sample='unif', injection_parameters=injection_parameters, 
    nlive=500, dlogz=4, clean=True, npool=16)


