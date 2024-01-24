import numpy as np

from astropy import constants 

def Planck(freq, temp): # freq: cm and temp: K 
# planck function in terms of frequency and temperature in cgs units 
    h = constants.h.cgs.value # cm^2 g s^-1
    c = constants.c.cgs.value # cm/s
    k = constants.k_B.cgs.value # cm^2 g s^-2 K^-1
    return (2*h*freq**3) / ((c**2)*(np.exp((h*freq)/(k*temp)) - 1)) # intensity: ergs/cm2/s/st/Hz

def box_integrator(x, y): # array of x and y points
# box rule integrator function 
    dx = x[2] - x[1]
    dy = y[1]
    area = dx*dy
    for i in range(2, len(x)-1):
        dx = x[i+1] - x[i]
        dy = y[i]
        area += dx*dy
    return area

def integrate_Planck(temp, x_min, x_max, density): # temp: K, wavenumber min and max value: μm^-1, density of wavenumber points
# function to integrate Planck 
    wavenumber_microns = np.arange(x_min, x_max, density)
    wavelengths = 1 / wavenumber_microns
    intensity_function = Planck(constants.c.cgs.value / (wavelengths/1e4), temp)
    return box_integrator(wavenumber_microns, intensity_function)

def intensity_integral(n, x):
# function for the integral portion of the emergent intensity equation 
    return (x**n)*(np.exp(-x))

def integrate_intensity(n, x_min, x_max, step):
# function to integrate the intensity integral
    x_array = np.arange(x_min, x_max, step)
    intensity_function = intensity_integral(n, x_array)
    return box_integrator(x_array, intensity_function)