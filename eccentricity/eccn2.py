from __future__ import division, print_function
import matplotlib.pyplot as plt
import bilby
import numpy as np

import hphc

C = 299792458.
G = 6.67408*1e-11
Mo = 1.989*1e30
Mpc = 3.086*1e22

sampling_frequency = 1024.0
duration = 4.0

def eccentric_waveform(frequency_array_, chirp_mass, mass_ratio, eccentricity, luminosity_distance, theta_jn, psi, phase, geocent_time, ra, dec):
    
    mass_1 = (chirp_mass*(1+mass_ratio)**(1/5))/mass_ratio**(3/5)
    mass_2 = chirp_mass*mass_ratio**(2/5)*(1+mass_ratio)**(1/5)
    luminosity_distance = luminosity_distance*Mpc
    total_mass = (mass_1+mass_2)*Mo
    symmetric_mass_ratio = (mass_1*mass_2)/((mass_1+mass_2)**2)
    maximum_f = (C**3)/( G*(mass_1+mass_2)*Mo*np.pi*6**(3/2) )
    minimum_f = 20.0
    mass_diff = (mass_1-mass_2)*Mo
    
    foo = np.array(frequency_array_)
    N = len(foo)
    h_plus = np.zeros(N)*1j
    h_cross = np.zeros(N)*1j

    mask_3 = np.logical_and(foo >= minimum_f, foo <= maximum_f)
    index = np.array(np.where(mask_3)).flatten()
    
    for k in index:
        arg_plus = {'iota_':theta_jn, 'beta_':psi, 'D_':luminosity_distance , \
                    'f_':foo[k], 'f0_':20.0, 'et0_':eccentricity, 'phic_':phase, \
                    'tc_':geocent_time, 'M_':total_mass, 'eta_':symmetric_mass_ratio, \
                    'ff_':maximum_f, 'delta_':mass_diff}

        fplus = hphc.Fn(**arg_plus)

        h_plus[k] = fplus.htilde()[0]
        h_cross[k] = fplus.htilde()[1]

    #return {'plus': h_plus, 'cross': h_cross}
    return {'plus': h_plus, 'cross': h_cross}

injection_parameters = dict(chirp_mass=13.488088173501811, mass_ratio=0.4, eccentricity=0.1, luminosity_distance=200.0, theta_jn=0.4, psi=0.1, phase=1.2, geocent_time=1180002601.0, ra=45.0, dec=5.73)

waveform_generator = bilby.gw.waveform_generator.WaveformGenerator(
    duration=duration, sampling_frequency=sampling_frequency,
    frequency_domain_source_model=eccentric_waveform)

minimum_frequency = 20.0+5.0
MM1 = 25.0
MM2 = 10.0
maximum_frequency = (C**3)/( G*(MM1+MM2)*Mo*np.pi*6**(3/2) )-5.0

ifos = bilby.gw.detector.InterferometerList(['H1', 'L1'])
for ifo in ifos:
    ifo.minimum_frequency = minimum_frequency
    ifo.maximum_frequency = maximum_frequency
ifos.set_strain_data_from_power_spectral_densities(
    sampling_frequency=sampling_frequency, duration=duration,
    start_time=injection_parameters['geocent_time'] - 3)
ifos.inject_signal(waveform_generator=waveform_generator,
                   parameters=injection_parameters)


prior = bilby.core.prior.PriorDict()
prior['chirp_mass'] = bilby.core.prior.Uniform(name='chirp_mass', minimum=5.0,maximum=30.0)
prior['mass_ratio'] = bilby.core.prior.Uniform(name='mass_ratio', minimum=0.1, maximum=1)
prior['eccentricity'] = bilby.core.prior.Uniform(name='eccentricity', minimum=0.01, maximum=0.2)

prior["luminosity_distance"] = 200.0
prior["theta_jn"] = 0.4
prior["psi"] = 0.1
prior["phase"] = 1.2
prior["geocent_time"] = 1180002601.0
prior["ra"] = 45.0
prior["dec"] = 5.73

'''
#this is to run the code faster
likelihood = bilby.gw.likelihood.GravitationalWaveTransient(
    interferometers=ifos, waveform_generator=waveform_generator, priors=prior,
    time_marginalization=False, phase_marginalization=False, distance_marginalization=False)

result_short = bilby.run_sampler(
    likelihood, prior, sampler='dynesty', outdir='short2', label="eccn2",
    nlive=1000, dlogz=3, npool=16,  # <- Arguments are used to make things fast - not recommended for general use
    clean=True
)
'''

likelihood = bilby.gw.likelihood.GravitationalWaveTransient(
    interferometers=ifos, waveform_generator=waveform_generator, priors=prior)

result_short = bilby.run_sampler(
    likelihood, prior, sampler='dynesty', outdir='short3', label="eccn3", npool=16, clean=True
)

























