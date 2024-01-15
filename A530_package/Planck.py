import numpy as np

from astropy import constants 

def Planck(freq, temp): # freq: cm and temp: K 
# planck function in terms of frequency and temperature in cgs units 
    h = constants.h.cgs.value # cm^2 g s^-1
    c = constants.c.cgs.value # cm/s
    k = constants.k_B.cgs.value # cm^2 g s^-2 K^-1
    return (2*h*freq**3) / ((c**2)*(np.exp((h*freq)/(k*temp)) - 1)) # intensity: ergs/cm2/s/st/Hz

def box_integrator(x, y):
    dx = x[2] - x[1]
    dy = y[2] - y[1]
    area = dx*dy
    for i in range(2, len(x)-1):
        dx = x[i+1] - x[i]
        dy = y[i+1] - y[i]
        area += dx*dy
    return area

def integrate(x, y, x_min, x_max):



"""
Then, wrap this integrator in another function that lets you specify the the minimum and
maximum x-axis value and density of points in the integration of a given function.
"""
