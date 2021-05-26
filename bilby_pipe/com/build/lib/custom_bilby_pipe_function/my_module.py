from bilby.gw.source import lal_binary_black_hole
from bilby.gw.waveform_generator import WaveformGenerator
import gwsurrogate

# Load the surrogate model
surrogate = gwsurrogate.LoadSurrogate('NRHybSur3dq8')

# Define the time-domain model
def my_function(times, mass, ratio, distance, theta, phi, speed, v_the, v_phi,**kwargs):

    arg_ = {'surrogate_':surrogate,'M_':mass, 'q_':ratio, 'dis_':distance, 'the_':theta, 'phi_':phi, 'v_mag_':speed, 'v_the_':v_the, 'v_phi_':v_phi, 'times_':times, **kwargs}
    
    fplus = WaVel.Fn(**arg_)
    
    h_plus = fplus.WaVe()[0]
    h_cross = fplus.WaVe()[1]

    return {'plus': h_plus, 'cross': h_cross}

def my_waveform_generator(duration, sampling_frequency, frequency_domain_source_model, **kwargs):
    
    return WaveformGenerator(duration=duration, sampling_frequency=sampling_frequency, frequency_domain_source_model=frequency_domain_source_model)

'''

# Python program to illustrate
# *args with first extra argument
def myFun(arg1, *argv):
    print ("First argument :", arg1)
    for arg in argv:
        print("Next argument through *argv :", arg)
 
myFun('Hello', 'Welcome', 'to', 'GeeksforGeeks')


# Python program to illustrate  **kargs for
# variable number of keyword arguments with
# one extra argument.
 
def myFun(arg1, **kwargs):
    for key, value in kwargs.items():
        print ("%s == %s" %(key, value))
 
# Driver code
myFun("Hi", first ='Geeks', mid ='for', last='Geeks') 

'''
