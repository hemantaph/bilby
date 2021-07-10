# written by Hemanta Phurailatpam (CUHK)
# python code for parameter estimation with eccentric waveform model (accuracy: 0PN at amplitude, 3PN et0^6 at fourier phase)

import bilby
import numpy as np

# calling my function from hphc.py file
import hphc

C = 299792458.
G = 6.67408*1e-11
Mo = 1.989*1e30
Mpc = 3.086*1e22

sampling_frequency = 1024.0
duration = 4.0

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
    #f_max is set according to the priors we choose for chirp_mass and mass_ratio; depends on max possible lso_f. 
    f_max = 110.
    #time of coalescence is taken to be 0 for convenience
    tc = 0.
    
    foo = np.array(frequency_array_, dtype='float')

    h_plus,h_cross = hphc.htilde(foo, total_mass, symmetric_mass_ratio, mass_diff, initial_eccentricity, luminosity_distance, theta_jn, psi, phase, tc, f_min, lso_f, f_max)

    return {'plus': h_plus, 'cross': h_cross}
#####################################################################

#'mass_1':23.5, 'mass_2':21.5
injection_parameters = dict(chirp_mass=19.564163812778446, mass_ratio=0.9148936170212766, initial_eccentricity=0.1, luminosity_distance=400.0, theta_jn=0.4, psi=0.1, phase=1.2, geocent_time=1180002601.0, ra=45.0, dec=5.73)

waveform_generator = bilby.gw.waveform_generator.WaveformGenerator(
    duration=duration, sampling_frequency=sampling_frequency,
    frequency_domain_source_model=eccentric_waveform)

minimum_frequency = 23.0
#f_max is set according to the priors we choose for chirp_mass and mass_ratio; depends on max possible lso_f. 
maximum_frequency = 107.0

ifos = bilby.gw.detector.InterferometerList(['H1', 'L1'])
for ifo in ifos:
    ifo.minimum_frequency = minimum_frequency
    ifo.maximum_frequency = maximum_frequency
ifos.set_strain_data_from_power_spectral_densities(
    sampling_frequency=sampling_frequency, duration=duration,
    start_time=injection_parameters['geocent_time'] - 3.5)
ifos.inject_signal(waveform_generator=waveform_generator,
                   parameters=injection_parameters)


prior = bilby.core.prior.PriorDict()
prior['chirp_mass'] = bilby.core.prior.Uniform(name='chirp_mass', minimum=17.411,maximum=21.764)
prior['mass_ratio'] = bilby.core.prior.Uniform(name='mass_ratio', minimum=0.5, maximum=1)
prior['initial_eccentricity'] = bilby.core.prior.LogUniform(name='eccentricity', minimum=0.01, maximum=0.2)
prior["luminosity_distance"] = 400.0
prior["theta_jn"] = 0.4
prior["psi"] = 0.1
prior["phase"] = 1.2
prior["geocent_time"] = 1180002601.0
prior["ra"] = 45.0
prior["dec"] = 5.73


#this is to run the code faster
likelihood = bilby.gw.likelihood.GravitationalWaveTransient(
    interferometers=ifos, waveform_generator=waveform_generator, priors=prior,
    time_marginalization=False, phase_marginalization=False, distance_marginalization=False)

result_short = bilby.run_sampler(
    likelihood, prior, sampler='dynesty', outdir='short', label="eccn",
    nlive=1000, dlogz=0.1, npool=8 )

























