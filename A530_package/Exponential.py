import numpy as np
from A530_package import box_integrator

def special_exp_integral(n, t):
# function for the special exponential integral portion 
    return (np.exp(-t))/(t**n)

def integrate_special_exp(n, x_min, x_max, step):
# function to integrate the special exponent 
    x_array = np.arange(x_min, x_max, step)
    integral = special_exp_integral(n, x_array)
    return (x_min**(n-1))*box_integrator(x_array, integral)