import numpy as np

#required constants
I = 13.598
loge = 0.43429
R = 1.09678*10**(-3)
alpha_0 = 1.0443*10**(-26)
a_0 = 0.1199654
a_1 = -1.18267*10**(-6)
a_2 = 2.64243*10**(-7)
a_3 = -4.40524*10**(-11)
a_4 = 3.23992*10**(-15)
a_5 = -1.39568*10**(-19)
a_6 = 2.78701*10**(-24)

def theta(T):
#returns theta value for a given temperature
    return 5040/T

def chi_n(n):
#returns chi_n for a given n
    return 13.598*(1-(1/n**2))

def chi_lambda(wavelength):
#returns chi_lambda for a given wavelength
    return (1.2398*10**4) / wavelength

def f_0(wavelength):
#returns f_0 for a given wavelength
    return -2.2763 - 1.685*np.log10(wavelength) + 0.76661*(np.log10(wavelength))**2 - 0.0533464*(np.log10(wavelength))**3

def f_1(wavelength):
#returns f_1 for a given wavelength
    return 15.2827 - 9.2846*np.log10(wavelength) + 1.99381*(np.log10(wavelength))**2 - 0.142631*(np.log10(wavelength))**3

def f_2(wavelength):
#returns f_2 for a given wavelength
    return -197.789 + 190.266*np.log10(wavelength) - 67.9775*(np.log10(wavelength))**2 + 10.6913*(np.log10(wavelength))**3 - 0.625151*(np.log10(wavelength))**4

def gaunt_ff(wavelength, temp):
#returns gaunt factor for ff for a given wavelength and temp
    return 1 + (0.3456/(wavelength*R)**(1/3))*((loge/(theta(temp)*chi_lambda(wavelength))) + 0.5)

def gaunt_bf(wavelength, n):
#returns gaunt factor for bf for a given wavelength and n 
    return 1 - ((0.3456/((wavelength*R)**(1/3)))*(((wavelength*R)/(n**2)) - 0.5))

def alpha_bf_neg_H(wavelength):
#returns alpha bf of negative H for a given wavelength
    return (a_0 + a_1*wavelength + a_2*wavelength**2 + a_3*wavelength**3 + a_4*wavelength**4 + a_5*wavelength**5 + a_6*wavelength**6)*(10**(-17))

def opacity_neg_H_bf(Pe, wavelength, temp):
#returns the opacity of negative H bf for a given Pe, wavelength, and temperature
    return 4.158*10**(-10)*alpha_bf_neg_H(wavelength)*Pe*(theta(temp)**(5/2))*(10**(0.754*theta(temp)))

def opacity_neg_H_ff(Pe, wavelength, temp):
#returns the opacity of negative H ff for a given Pe, wavelength, and temperature
    return 10**(-26)*Pe*(10**(f_0(wavelength) + f_1(wavelength)*np.log10(theta(temp)) + f_2(wavelength)*(np.log10(theta(temp)))**2))

def opacity_H_ff(wavelength, temp):
#returns the opacity of neutral H ff for a given wavelenth and temperature
    return alpha_0*(wavelength**3)*gaunt_ff(wavelength,temp)*(loge/(2*theta(temp)*I))*(10**(-theta(temp)*I))

def opacity_H_bf(wavelength, temp):
#returns the opacity of neutral H bf for a given wavelenth and temperature
    total = 0
    for n in range(1, 100):
        if wavelength > (n**2/R):
            total = total + 0
        else:
            total = total + ((wavelength/n)**3)*gaunt_bf(wavelength,n)*(10**(-theta(temp)*chi_n(n)))
    return total*alpha_0

def opacity_total(H_bf, H_ff, neg_H_bf, neg_H_ff, temp, wavelength, Pe):
#returns total opacity for a given Pe, wavelength, and temperature
    return ((H_bf + H_ff + neg_H_bf)*(1-10**(-theta(temp)*chi_lambda(wavelength))) + neg_H_ff)
