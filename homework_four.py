from A530_package import quadratic_source, linear_source, eddington_flux
from scipy import integrate
import scipy.special as sc
import numpy as np

"""
7. Emergent Eddington Flux Hv(0)
Use your work from Problem 5 to write a function that calculates Hv(0) given an array of
optical depths τv using an exponential integral. Choose a quadratic form for Sv(τv) and compute Hv (0) for it.
Confirm that this works as you expect it to for a linear source function.
"""

#all related functions are found under A530_package/Exponential.py

#hardwire values for a_n coefficients
a0 = a1 = 1
a2 = 2

#set array of optical depths 
x_min = 10**(-5) 
x_max = 30
step = 10**(-6) 
tau_array = np.power(10, np.arange(np.log10(x_min), np.log10(x_max), step))

print("eddington flux for quadratic source function: {}".format(eddington_flux(quadratic_source, tau_array, a0, a1, a2)))
#eddington flux for quadratic source function: 0.9166616669497857
print("eddington flux for linear source function: {}".format(eddington_flux(linear_source, tau_array, a0, a1, a2)))
#eddington flux for linear source function: 0.41666166695251616

#smaller eddington flux for linear source function as expected 