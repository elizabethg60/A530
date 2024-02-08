import numpy as np
import scipy.special as sc
from scipy import integrate
from A530_package import box_integrator

def special_exp_integral(n, t):
# function for the special exponential integral portion 
    return (np.exp(-t))/(t**n)

def integrate_special_exp(n, x_min, x_max, step):
# function to integrate the special exponent 
    x_array = np.arange(x_min, x_max, step)
    integral = special_exp_integral(n, x_array)
    return (x_min**(n-1))*box_integrator(x_array, integral)

def quadratic_source(tau_array, a0, a1, a2):
# returns a quadratic source function for given a_n values and array of optical depth
    return a0 + a1*tau_array + a2*tau_array**2

def linear_source(tau_array, a0, a1, a2):
# returns a linear source function for given a_n values and array of optical depth
    return a0 + a1*tau_array

def eddington_flux(source_func, tau_array, a0, a1, a2):
# computes the eddington flux 
    E_2 = sc.expn(2, tau_array)
    S_v = source_func(tau_array, a0, a1, a2)
    return 0.5*integrate.simpson(E_2*S_v, tau_array)

def eddington_flux_planck(planck, tau_array):
# computes the eddington flux under LTE when source = planck 
    E_2 = sc.expn(2, tau_array)
    S_v = planck
    return 0.5*integrate.simpson(E_2*S_v, tau_array)

def D(func, x, h):
# computes the numerical derivative 
    return (func(x+h)-func(x))/h