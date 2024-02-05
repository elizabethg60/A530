from A530_package import Planck, eddington_flux_planck
import matplotlib.pyplot as plt
from astropy import constants
import numpy as np

"""
(b) Modify your function from Problem 7 to calculate Fν (0) as a function of wavenumber
for a Teff =8,700 K grey atmosphere. Use a wavenumber range similar to that in Problem
2.
"""
#set array of optical depths 
x_min = 10**(-16) 
x_max = 50
step = 10**(-5) 
tau_array = np.power(10, np.arange(np.log10(x_min), np.log10(x_max), step))

#set array of wavenumber 
wav = np.logspace(np.log10(1000), np.log10(110000), num = len(tau_array))

Teff = 8700
#temperature in terms of Teff and tau (part a)
T_fct_tau = Teff * ((3/4)*(tau_array+(2/3)))**(1/4)

#planck function in terms of wavenumber and tau
planck_function = Planck(constants.c.cgs.value / (wav/1e8), T_fct_tau)

print("surface flux is {}".format(4*np.pi*eddington_flux_planck(planck_function, tau_array)))
