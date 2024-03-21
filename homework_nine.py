import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from A530_package import Pe_converge, saha

"""
18. The VALiiiC Atmosphere
a) Write a function that interprets VALIIIC.txt for use with your other programs.
"""

#read in data
path = '/Users/efg5335/Desktop/Courses/A530/data/'
VALIIIC = pd.read_csv(path + 'Solar_Data/VALIIIC_sci_e.txt', sep= ' ', header=None)
VALIIIC.set_axis(["h", "m", "tau_500", "T", "V", "n_H", "n_e", "Ptotal", "Pgas/Ptotal", "rho"], axis=1, inplace=True)
# h     m     tau_500     T     V     n_H     n_e     Ptotal  Pgas/Ptotal  rho
#(km) (g cm )            (K)  (km/s) (cm^-3) (cm^-3) (dyn cm^-2)          (g cm^-3)

"""
b) Using your Φ function for hydrogen and the electron gas pressure implied by ne and
T in the VALiiiC atmosphere, reproduce the top part of Figure 8.8 (which comes from
the VALIII paper). 

Do you get the right number of protons? Does your answer at least
make physical sense compared to the number of electrons? Remember that Vernazza et
al. used NLTE methods to generate that figure — what ionization temperature do you
need to assume at h = 800km to reproduce their np there?
"""

temperature = VALIIIC["T"]
Pg = VALIIIC["Ptotal"]*VALIIIC["Pgas/Ptotal"]
n_H = VALIIIC["n_H"]

Pe_array = []
n_p = []
for i in range(0, len(temperature)):
    #initial guess for Pe 
    partition_function = pd.read_csv(path + 'RepairedPartitionFunctions.txt', sep= ' ', header=None)
    ionization_potential = pd.read_fwf(path + 'ioniz.txt', header=None)
    element_ion_pot = float(ionization_potential.loc[ionization_potential[1] == "H"][3])
    theta_arr = np.arange(0.2, 2.2, step = 0.2)
    temp_arr = 5040/theta_arr
    phi = saha("H", temperature[i], partition_function, temp_arr, element_ion_pot)
    Pe_initial = np.sqrt(phi*Pg[i])

    #converge to a Pe value 
    Pe_array.append(Pe_converge(Pg[i], temperature[i], Pe_initial))
    #determine corresponding n_p value
    n_p.append(n_H[i]*(phi/Pe_array[i]))

#reproduce the top part of Figure 8.8
plt.plot(VALIIIC["h"], np.log10(VALIIIC["n_e"]), label = '$n_{e}$', color = 'k')
plt.plot(VALIIIC["h"], np.log10(n_p), label = '$n_{p}$', color = 'r')
plt.xlim(0,800)
plt.ylim(8,14)
plt.gca().invert_xaxis()
plt.xlabel("h (km)")
plt.ylabel("log n")
plt.legend()
plt.savefig("Figures/hw_nine_figures/fig8_8.pdf")
plt.show()