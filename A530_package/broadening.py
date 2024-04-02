import numpy as np
import scipy

#required constants (cgs)
c = 3 * 10**10
h = 6.62 * 10**-27
k = 1.38 * 10**-16
m = 3.82 * 10**-23
micro = 0 

def B_ul(Aul, freq_center):
#returns Bul value given relation to Aul
    return Aul / ((2*h*freq_center**3)/ c**2)

def B_lu(Aul, freq_center, gu, gl):
#returns Blu value given statistical weight ratio from Bul
    Bul = B_ul(Aul, freq_center)
    return Bul * (gu/gl)

def doppler_width(freq_center, temp):
#returns doppler width
    return (freq_center/c)*np.sqrt(((2*k*temp)/m) + micro**2)

def C6_fct(I, freq_center):
#returns C6 constant (Gray 11.30) #chi_n = 0
    return 0.3*10**-30 * ((I-(h*freq_center*6.242*10**11))**-2 - (I)**-2)

def gamma_6(I, freq_center, Pg, temp):
#returns van der waals broadening (Gray 11.29)
    C6 = C6_fct(I, freq_center)
    return 10**(20 + (2/5)*np.log10(C6) + np.log10(Pg) - (7/10)*np.log10(temp))

def gamma_4(C4, Pe, temp):
#returns quadratic stark broadening (Gray 11.27)
    return 10**(19 + (2/3)*C4 + np.log10(Pe) - (5/6)*np.log10(temp))

def gamma_total(Aul, C4, Pe, temp, I, freq_center, Pg):
    stark = gamma_4(C4, Pe, temp)
    van = gamma_6(I, freq_center, Pg, temp)
    return 4*np.pi*Aul + stark + van 

def u(freq, freq_center, dopplerwidth):
#returns u value for the Voigt profile
    return (freq - freq_center) / dopplerwidth

def a(gamma, dopplerwidth):
#returns a value for the Voigt profile
    return gamma / (4*np.pi*dopplerwidth)

def extinction_coefficient(Aul, gu, gl, freq_center, freq, temp, I, C4, Pe, Pg):
#returns the extinction coefficient
    dopplerwidth = doppler_width(freq_center, temp)
    gamma = gamma_total(Aul, C4, Pe, temp, I, freq_center, Pg)

    a_constant = a(gamma, dopplerwidth)
    u_constant = u(freq, freq_center, dopplerwidth)
    voigt = np.real(scipy.special.wofz(u_constant + 1j * a_constant))

    Blu = B_lu(Aul, freq_center, gu, gl)

    return ((h*freq_center)/(4*np.pi)) * Blu * voigt