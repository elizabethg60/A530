import math
import numpy as np
import matplotlib.pyplot as plt

from A530_package import Planck
from astropy import units as u
from astropy import constants

"""
(a) Write your own Planck function you trust that accepts vector input and in units you
understand. Use this to confirm and understand the input and output of a library Planck
function. In future assignments, you may use your own function or a library function.
"""

from astropy.modeling.models import BlackBody
bb = BlackBody(temperature=5778*u.K)
wav = np.arange(1000, 110000) * u.AA
intensity_astropy = bb(wav)
plt.plot(wav, intensity_astropy, label = "astropy")

intensity_function = Planck(constants.c.cgs.value / (wav.value/1e8), 5778)
plt.scatter(wav, intensity_function, label = "my function", s = 1, color = 'r')
plt.legend()
plt.savefig("Figures/hw_one_figures/sample_a.png")
plt.show()

"""
b) Make three plots: one of Bν vs.  ̃ν on the range 0 to 12 μm−1, one of log Bν vs.  ̃ν on
the same range, and one of log Bν vs. log  ̃ν (on the range −1.0 < log  ̃ν < 1.2).

On all 3 plots, show curves for T=10,000 K, 7,000 K, and 3,000 K. “log” is base 10.
Choose appropriately illustrative y ranges, and include a legend
"""

#wavenumber range given
wavenumber_microns = np.arange(0, 12, .2)
#wavelength range 
wavelengths = 1 / wavenumber_microns
intensity_function3000 = Planck(constants.c.cgs.value / (wavelengths/1e4), 3000)
intensity_function7000 = Planck(constants.c.cgs.value / (wavelengths/1e4), 7000)
intensity_function10000 = Planck(constants.c.cgs.value / (wavelengths/1e4), 10000)
plt.plot(wavenumber_microns, intensity_function3000, label = "3000K")
plt.plot(wavenumber_microns, intensity_function7000, label = "7000K")
plt.plot(wavenumber_microns, intensity_function10000, label = "10000K")
plt.legend()
plt.savefig("Figures/hw_one_figures/sample_bi.png")
plt.show()

plt.plot(wavenumber_microns, intensity_function3000, label = "3000K")
plt.plot(wavenumber_microns, intensity_function7000, label = "7000K")
plt.plot(wavenumber_microns, intensity_function10000, label = "10000K")
plt.legend()
plt.yscale("log")
plt.savefig("Figures/hw_one_figures/sample_bii.png")
plt.show()

#wavenumber range given
wavenumber_microns = np.linspace(10**(-1.0), 10**(1.2), 50)
#wavelength range 
wavelengths = 1 / wavenumber_microns
intensity_function3000 = Planck(constants.c.cgs.value / (wavelengths/1e4), 3000)
intensity_function7000 = Planck(constants.c.cgs.value / (wavelengths/1e4), 7000)
intensity_function10000 = Planck(constants.c.cgs.value / (wavelengths/1e4), 10000)
plt.plot(wavenumber_microns, intensity_function3000, label = "3000K")
plt.plot(wavenumber_microns, intensity_function7000, label = "7000K")
plt.plot(wavenumber_microns, intensity_function10000, label = "10000K")
plt.legend()
plt.yscale("log")
plt.xscale("log")
plt.savefig("Figures/hw_one_figures/sample_biii.png")
plt.show()

#make plots pretty
#finalize code (last read through of homework) and save everything
#read latex guidelines and write up latex portion
#submit! 