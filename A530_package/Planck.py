import numpy as np

from astropy import constants 

def Planck(freq, temp): # freq: cm and temp: K 
# planck function in terms of frequency and temperature in cgs units 
    h = constants.h.cgs.value # cm^2 g s^-1
    c = constants.c.cgs.value # cm/s
    k = constants.k_B.cgs.value # cm^2 g s^-2 K^-1
    return (2*h*freq**3) / ((c**2)*(np.exp((h*freq)/(k*temp)) - 1)) # intensity: ergs/cm2/s/st/Hz
