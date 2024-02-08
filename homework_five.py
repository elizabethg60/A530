from A530_package import Planck, eddington_flux_planck, D
import matplotlib.pyplot as plt
from astropy import constants
import numpy as np

"""
(b) Modify your function from Problem 7 to calculate Fν (0) as a function of wavenumber
for a Teff =8,700 K grey atmosphere. Use a wavenumber range similar to that in Problem 1.
(c) What is the Eddington-Barbier expectation for Fν (0) given this temperature distri-
bution? By this I mean: the proper calculation to determine Fν (0) from T (τ ) involves a
non-analytic integral, but you can approximate the answer using the Eddington-Barbier
expression for Fν (0). 
"""

#set array of optical depths 
x_min = 10**(-4) 
x_max = 50
step = 10**(-5) 
tau_array = np.power(10, np.arange(np.log10(x_min), np.log10(x_max), step))

#wavenumber range given in hw 1
wavenumber_microns = np.arange(0, 12, .05)
#wavelength range 
wavelengths = 1 / wavenumber_microns

Teff = 8700
#temperature in terms of Teff and tau (part a)
T_fct_tau = Teff * ((3/4)*(tau_array+(2/3)))**(1/4)

planck_function = []
EB_approx = []
for wav in wavelengths:
    #planck function in terms of wavenumber and tau
    planck_function.append(Planck(constants.c.cgs.value / (wav/1e4), T_fct_tau))
    #EB approx for tau = 2/3 
    EB_approx.append(Planck(constants.c.cgs.value / (wav/1e4), Teff * ((3/4)*((2/3)+(2/3)))**(1/4))*np.pi)

#get emergent flux from planck function 
flux_arr = 4*np.pi*eddington_flux_planck(planck_function, tau_array)

"""
(d) Plot Fν (0) vs.  ̃ν, comparing your the answers in parts (b) and (c). Where does the
E-B approximation underestimate and overestimate the true emergent flux? What does
this tell you about the functional form of the source function at these frequencies? (Hint:
consider problem 2f.)
"""

plt.plot(wavenumber_microns, np.array(flux_arr)/(10**(-4)), linewidth = 2, label = 'true', color = 'r')
plt.plot(wavenumber_microns, np.array(EB_approx)/(10**(-4)), linewidth = 2, label = 'EB', color = 'k')
plt.xlabel("wavenumber [$μm^{−1}$]", fontsize = 12)
plt.ylabel("emergent flux [$10^{-4}$ erg $cm^{-2}$ $s^{-1}$ $Hz^{-1}$]", fontsize = 12)
plt.xticks(fontsize = 12) 
plt.yticks(fontsize = 12) 
plt.legend(fontsize = 12)
plt.legend()
plt.savefig("Figures/hw_five_figures/sample_d.pdf", bbox_inches='tight')
plt.show()

"""
(e) Numerically calculate d2Sν /dτ 2 at τ = 1 (this is your proxy for finding the the sign of
a2; “numerically calculate” means take finite differences to approximate the derivative).
Use these values as a function of  ̃ν to quantitatively illustrate your answer to (d) (that is,
plot the difference between the exact and E-B answers as a function of wavenumber, and
show that this difference is what you expect, given your estimate of the sign of a2 from
finite differencing). 
"""

def Planck_tau(tau):
# Planck function in terms of tau
    return (Planck(constants.c.cgs.value / (wav/1e4), Teff * ((3/4)*((tau)+(2/3)))**(1/4))*np.pi)

h = 10**(-6) 
second_D = []
for wav in wavelengths:
    prime_one = D(Planck_tau, 1, h)
    prime_oneh = D(Planck_tau, 1+h, h)
    second_D.append((prime_oneh-prime_one)/h)

plt.plot(wavenumber_microns, flux_arr - EB_approx, linewidth = 2, label = 'exact - EB', color = 'b')
plt.scatter(wavenumber_microns, second_D, linewidth = 2, label = 'second derivative', color = 'r')
plt.xlabel("wavenumber [$μm^{−1}$]", fontsize = 12)
plt.xticks(fontsize = 12) 
plt.yticks(fontsize = 12) 
plt.legend(fontsize = 12)
plt.legend()
plt.savefig("Figures/hw_five_figures/sample_e.pdf", bbox_inches='tight')
plt.show()