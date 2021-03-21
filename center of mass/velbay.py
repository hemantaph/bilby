""" Bayes analysis for moving sources """

import sys
import subprocess

import bilby
import numpy as np


# Define labels, sampling frequency and times
outdir = 'outdir'
label = 'moving_bbh'
sampling_frequency = 1024
duration = 2.0
geocent_time = 0.0
start_time = geocent_time - duration + 0.02


# Define the time-domain model
def moving_bbh(times, mass, ratio, distance, theta, phi, speed, v_the, v_phi):
	# Transform times into a string
	times_str = ' '.join(map(str,times))

	# Run 'WaVel.py' and import the time and polarizations
	subprocess.run([sys.executable, "WaVel.py", str(mass), str(ratio), str(distance), str(theta), str(phi), str(speed), str(v_the), str(v_phi), times_str])

	# 'Get' the polarizations
	WaVel = np.load('WaVel.npy')


	# Return results
	return {'plus': WaVel[0], 'cross': WaVel[1]}
	

# Define parameters to inject
injection_parameters = dict(mass=50, ratio=6, distance=300, theta=45, phi=45, speed=3000, v_the=0, v_phi=0, ra=0, dec=0, psi=0, geocent_time=geocent_time)


# Generate model with waveform_generator
waveform = bilby.gw.waveform_generator.WaveformGenerator(sampling_frequency=sampling_frequency, duration=duration, time_domain_source_model=moving_bbh, start_time=start_time)  


# Inject signal into ifos
ifos = bilby.gw.detector.InterferometerList(['H1', 'L1', 'V1'])
ifos.set_strain_data_from_power_spectral_densities(sampling_frequency=sampling_frequency, duration=duration, start_time=start_time)
ifos.inject_signal(waveform_generator=waveform, parameters=injection_parameters)


# Create priors
prior = injection_parameters.copy()
prior['mass'] = bilby.core.prior.Uniform(45, 55, r'$M$', unit='$M_sun$')
#prior['ratio'] = bilby.core.prior.Uniform(5, 7, r'$q$')
#prior['distance'] = bilby.core.prior.Uniform(275, 325, r'd', unit='$Mpc$')
#prior['theta'] = bilby.core.prior.Uniform(40, 50, r'$\theta$', unit='degree')
#prior['phi'] = bilby.core.prior.Uniform(40, 50, r'$\phi$', unit='degree')
#prior['speed'] = bilby.core.prior.Uniform(0, 3000, r'$|v|$', unit='$km/s$')
#prior['v_the'] = bilby.core.prior.Uniform(0, 40, r'$v_\theta$', unit='$degree$')
#prior['v_phi'] = bilby.core.prior.Uniform(0, 40, r'$v_\phi$', unit='$degree$')


# Define likelihood
likelihood = bilby.gw.likelihood.GravitationalWaveTransient(interferometers=ifos, waveform_generator=waveform, priors=prior)


# Launch sampler
result = bilby.core.sampler.run_sampler(likelihood, prior, sampler='dynesty', npoints=200, injection_parameters=injection_parameters, outdir=outdir, label=label, dlogz=3, npool=6)


# Plot results
result.plot_corner()
