from A530_package import Planck, box_integrator
import matplotlib.pyplot as plt
from astropy import constants
import numpy as np

#set temperature (Kelvin)
temp = 7500
#set wavenumber range: 0 to "inf"
wavenumber_microns = np.linspace(0, 1e5, int(1e7))
#wavelength range 
wavelengths = 1 / wavenumber_microns
intensity_function = Planck(constants.c.cgs.value / (wavelengths/1e4), temp)
plt.plot(wavenumber_microns, intensity_function * 1e4, label = "{}K".format(temp), linewidth = 5)
plt.xlabel("wavenumber [$μm^{−1}$]", fontsize = 12)
plt.ylabel("intensity [$10^{-4}$ erg $cm^{-2}$ $s^{-1}$ $st^{-1}$ $Hz^{-1}$]", fontsize = 12)
plt.xticks(fontsize = 12) 
plt.yticks(fontsize = 12) 
plt.legend(fontsize = 12)
plt.xscale("log")
plt.show()

print("area under curve: {}".format(box_integrator(wavenumber_microns, intensity_function * 1e4))) #= -2.0511731667056087e-06 (units of 1e-4 specific intensity times wavenumber)

"""
Then, wrap this integrator in another function that lets you specify the the minimum and
maximum x-axis value and density of points in the integration of a given function.
"""