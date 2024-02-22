from astropy import constants
import numpy as np

#constants
alpha_0 = 1.0443*10**(-26)
I = 12.598
a_0 = 0.1199654
a_1 = -1.18267*10**(-6)
a_2 = 2.64243*10**(-7)
a_3 = -4.40524*10**(-11)
a_4 = 3.23992*10**(-15)
a_5 = -1.39568*10**(-19)
a_6 = 2.78701*10**(-24)

def chi_n(n):
    return 13.598*(1-(1/n**2))

def theta(T):
    return 5040/T

def alpha_cb_neg_H(wavelength):
    return a_0 + a_1*wavelength + a_2*wavelength**2 + a_3*wavelength**3 + a_4*wavelength**4 + a_5*wavelength**5 + a_6*wavelength**6

def f_0(wavelength):
    return -2.2763 - 1.685*np.log10(wavelength) + 0.76661*(np.log10(wavelength))**2 - 0.0533464*(np.log10(wavelength))**3

def f_1(wavelength):
    return 15.2827 - 9.2846*np.log10(wavelength) + 1.99381*(np.log10(wavelength))**2 - 0.142631*(np.log10(wavelength))**3

def f_2(wavelength):
    return -197.789 + 190.266*np.log10(wavelength) - 67.9775*(np.log10(wavelength))**2 + 10.6913*(np.log10(wavelength))**3 - 0.625151*(np.log10(wavelength))**4

