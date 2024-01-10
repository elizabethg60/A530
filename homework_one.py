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
plt.plot(wav, intensity_astropy/1e-5, label = "astropy", linewidth = 5)

intensity_function = Planck(constants.c.cgs.value / (wav.value/1e8), 5778)
plt.plot(wav, intensity_function/1e-5, label = "my function", color = 'r')

plt.xlabel("wavelength [Å]", fontsize = 12)
plt.ylabel("intensity [$10^{-5}$ erg $cm^{-2}$ $s^{-1}$ $st^{-1}$ $Hz^{-1}$]", fontsize = 12)
plt.xticks(fontsize = 12) 
plt.yticks(fontsize = 12) 
plt.legend(fontsize = 12)
plt.savefig("Figures/hw_one_figures/sample_a.png", bbox_inches='tight')
plt.show()

import pandas as pd
df = pd.DataFrame()
df["Wavelength"] = wav.value[0:10]
df["Astropy"] = intensity_astropy.value[0:10]
df["My Function"] = intensity_function[0:10]
with open('table.tex', 'w') as f:
    f.write(df.to_latex(index = False))

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
plt.plot(wavenumber_microns, intensity_function3000 *1e4, label = "3000K", linewidth = 5)
plt.plot(wavenumber_microns, intensity_function7000*1e4, label = "7000K", linewidth = 5)
plt.plot(wavenumber_microns, intensity_function10000*1e4, label = "10000K", linewidth = 5)
plt.xlabel("wavenumber [$μm^{−1}$]", fontsize = 12)
plt.ylabel("intensity [$10^{-4}$ erg $cm^{-2}$ $s^{-1}$ $st^{-1}$ $Hz^{-1}$]", fontsize = 12)
plt.xticks(fontsize = 12) 
plt.yticks(fontsize = 12) 
plt.legend(fontsize = 12)
plt.savefig("Figures/hw_one_figures/sample_bi.png", bbox_inches='tight')
plt.show()

plt.plot(wavenumber_microns, intensity_function3000, label = "3000K", linewidth = 5)
plt.plot(wavenumber_microns, intensity_function7000, label = "7000K", linewidth = 5)
plt.plot(wavenumber_microns, intensity_function10000, label = "10000K", linewidth = 5)
plt.xlabel("wavenumber [$μm^{−1}$]", fontsize = 12)
plt.ylabel("log intensity [erg $cm^{-2}$ $s^{-1}$ $st^{-1}$ $Hz^{-1}$]", fontsize = 12)
plt.xticks(fontsize = 12) 
plt.yticks(fontsize = 12) 
plt.legend(fontsize = 12)
plt.yscale("log")
plt.savefig("Figures/hw_one_figures/sample_bii.png", bbox_inches='tight')
plt.show()

#wavenumber range given
wavenumber_microns = np.linspace(10**(-1.0), 10**(1.2), 50)
#wavelength range 
wavelengths = 1 / wavenumber_microns
intensity_function3000 = Planck(constants.c.cgs.value / (wavelengths/1e4), 3000)
intensity_function7000 = Planck(constants.c.cgs.value / (wavelengths/1e4), 7000)
intensity_function10000 = Planck(constants.c.cgs.value / (wavelengths/1e4), 10000)
plt.plot(wavenumber_microns, intensity_function3000, label = "3000K", linewidth = 5)
plt.plot(wavenumber_microns, intensity_function7000, label = "7000K", linewidth = 5)
plt.plot(wavenumber_microns, intensity_function10000, label = "10000K", linewidth = 5)
plt.xlabel("log wavenumber [$μm^{−1}$]", fontsize = 12)
plt.ylabel("log intensity [erg $cm^{-2}$ $s^{-1}$ $st^{-1}$ $Hz^{-1}$]", fontsize = 12)
plt.xticks(fontsize = 12) 
plt.yticks(fontsize = 12) 
plt.legend(fontsize = 12)
plt.yscale("log")
plt.xscale("log")
plt.savefig("Figures/hw_one_figures/sample_biii.png", bbox_inches='tight')
plt.show()


#read latex guidelines and write up latex portion
#submit! 