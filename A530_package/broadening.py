import numpy as np
import scipy

#required constants (cgs)
first_lambda = 5889.95
first_g_lower = 2
first_g_upper = 4
first_A = 6.16 * 10**-1 * 10**8
first_stark = np.log10(-15.17)

second_lambda = 5895.924 
second_g_lower = 2
second_g_upper = 2
second_A = 6.14 * 10**-1 * 10**8
second_stark = np.log10(-15.33)

c = 3 * 10**10
h = 6.62 * 10**-27
k = 1.38 * 10**-16
m = 3.82 * 10**-23
micro = #TO DO 

def B_ul(Aul, freq):
#returns Bul value given relation to Aul
    return Aul / ((2*h*freq**3)/ c**2)

def B_lu(Bul, gu, gl):
#returns Blu value given statistical weight ratio from Bul
    return Bul * (gu/gl)

def doppler_width(freq, temp):
#returns doppler width
    return (freq/c)*np.sqrt((2*k*temp)/m + micro**2)

def v(freq_vary, freq, dopplerwidth):
#returns v value for the Voigt profile
    return (freq_vary - freq) / dopplerwidth

def a(gamma, dopplerwidth):
#returns a value for the Voigt profile
    return gamma / (4*np.pi*dopplerwidth)

def C6_fct(I, chi_n, freq):
#returns C6 constant (Gray 11.30)
    return 0.3*10**-30 * ((I-chi_n-h*freq)**-2 - (I - chi_n)**-2)

def gamma_4(C4, Pe, T):
#returns quadratic stark broadening (Gray 11.27)
    return np.log10(19 + (2/3)*np.log10(C4) + np.log10(Pe) - (5/6)*np.log10(T))

def gamma_6(C6, Pg, T):
#returns van der waals broadening (Gray 11.29)
    return np.log10(20 + (2/5)*np.log10(C6) + np.log10(Pg) - (7/10)*np.log10(T))

def gamma_total(natural, stark, van):
    return natural + stark + van 

def extinction_coefficient(freq, Blu, a, u):
#returns the extinction coefficient
    voigt = numpy.real(scipy.special.wofz(u + 1j * a))

    return ((h*freq)/(4*np.pi)) * Blu * voigt