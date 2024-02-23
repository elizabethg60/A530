from A530_package import opacity_total, opacity_neg_H_bf, opacity_neg_H_ff, opacity_H_ff, opacity_H_bf, saha
import matplotlib.pyplot as plt
import numpy as np

"""
a) Compute the total continuous opacity from hydrogen at a given electron pressure, temperature, and wavelength. 
b) As a check that you've got all of the opacity sources correct, reproduce Fig. 8.5(b-d).
"""

#hardwire temperature and electron pressure according to figure in Gray
temp = 5143
Pe = 10**1.08

#collect the values for each opacity term 
opacity_neg_H_bf_arr = []
opacity_neg_H_ff_arr = []
opacity_H_ff_arr = []
opacity_H_bf_arr = []

#iterate through wavelengths of 4000 - 20000
wavelength_arr = []
wavelength = 4000
neg = False
while wavelength < 20000:
    wavelength_arr.append(wavelength)

    neg_H_bf_value = opacity_neg_H_bf(Pe, wavelength, temp)
    if neg_H_bf_value >= 0.0 and neg == False:
        opacity_neg_H_bf_arr.append(opacity_neg_H_bf(Pe, wavelength, temp)/Pe)
    else: 
        opacity_neg_H_bf_arr.append(0.0)
        neg = True

    opacity_neg_H_ff_arr.append(opacity_neg_H_ff(Pe, wavelength, temp)/Pe)
    opacity_H_ff_arr.append(opacity_H_ff(wavelength, temp)/Pe)
    opacity_H_bf_arr.append(opacity_H_bf(wavelength, temp)/Pe)

    wavelength = wavelength + 1

#compute total continuous opacity from hydrogen given electron pressue, temperature, and wavelength
opacity_total_arr = opacity_total(np.array(opacity_H_bf_arr), np.array(opacity_H_ff_arr), np.array(opacity_neg_H_bf_arr), np.array(opacity_neg_H_ff_arr), temp, wavelength, Pe)

#reproduce Gray's figures
plt.plot(wavelength_arr,opacity_total_arr/10**(-26), label = "Total", color = 'k', linewidth = 5)
plt.plot(wavelength_arr,np.array(opacity_neg_H_bf_arr)/10**(-26), label = "$H^{-}_{bf}$", color = 'b', linewidth = 5, linestyle ="--")
plt.plot(wavelength_arr,np.array(opacity_neg_H_ff_arr)/10**(-26), label = "$H^{-}_{ff}$", color = 'r', linewidth = 5, linestyle ="-.")
plt.plot(wavelength_arr,np.array(opacity_H_ff_arr)/10**(-26), label = "$H_{ff}$", color = 'g', linewidth = 5, linestyle =":")
plt.plot(wavelength_arr,np.array(opacity_H_bf_arr)/10**(-26), label = "$H_{bf}$", color = 'y', linewidth = 5, linestyle ="dotted")
plt.xlabel("wavelength [Å]", fontsize = 12)
plt.ylabel("opacity / Pe [$10^{-26}$ $cm^{2}$ / H atom per dyne/$cm^{2}$", fontsize = 12)
plt.xticks(fontsize = 12) 
plt.yticks(fontsize = 12) 
plt.legend(fontsize = 12)
plt.savefig("Figures/hw_seven_figures/sample_bi.pdf", bbox_inches='tight')
plt.show()