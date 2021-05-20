import WaVel


def my_function(frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
                phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, **kwargs):
    
   # Load the surrogate model
surrogate = gwsurrogate.LoadSurrogate('NRHybSur3dq8')

# Define the time-domain model
def moving_bbh(times, mass, ratio, distance, theta, phi, speed, v_the, v_phi):

    arg_ = {'surrogate_':surrogate,'M_':mass, 'q_':ratio, 'dis_':distance, 'the_':theta, 'phi_':phi, 'v_mag_':speed, 'v_the_':v_the, 'v_phi_':v_phi, 'times_':times}
    
    h_plus,h_cross = WaVel.WaVe(**arg_)

    return {'plus': h_plus, 'cross': h_cross}
