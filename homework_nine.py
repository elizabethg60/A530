import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from A530_package import Pe_converge, saha

"""
18. The VALiiiC Atmosphere
a) Write a function that interprets VALIIIC.txt for use with your other programs.
"""

#read in data
path = '/Users/efg5335/Desktop/Courses/A530/data/'
solar_abundance = pd.read_csv(path + 'Solar_Data/SolarAbundance.txt', sep="	")
partition_function = pd.read_csv(path + 'RepairedPartitionFunctions.txt', sep= ' ', header=None)
ionization_potential = pd.read_fwf(path + 'ioniz.txt', header=None)
theta_arr = np.arange(0.2, 2.2, step = 0.2)
temp_arr = 5040/theta_arr

VALIIIC = pd.read_csv(path + 'Solar_Data/VALIIIC_sci_e.txt', sep= ' ', header=None)
VALIIIC.set_axis(["h", "m", "tau_500", "T", "V", "n_H", "n_e", "Ptotal", "Pgas/Ptotal", "rho"], axis=1, inplace=True)
# h     m     tau_500     T     V     n_H     n_e     Ptotal  Pgas/Ptotal  rho
#(km) (g cm )            (K)  (km/s) (cm^-3) (cm^-3) (dyn cm^-2)          (g cm^-3)

"""
b) Using your Φ function for hydrogen and the electron gas pressure implied by ne and
T in the VALiiiC atmosphere, reproduce the top part of Figure 8.8 (which comes from
the VALIII paper). Do you get the right number of protons? Does your answer at least
make physical sense compared to the number of electrons? Remember that Vernazza et
al. used NLTE methods to generate that figure — what ionization temperature do you
need to assume at h = 800km to reproduce their np there?
"""

temperature = VALIIIC["T"]
Pg = VALIIIC["Ptotal"]*VALIIIC["Pgas/Ptotal"]
n_H = VALIIIC["n_H"]
n_e = VALIIIC["n_e"]

# Pe_array = []
# n_p = []
# for i in range(0, len(temperature)):
#     #initial guess for Pe 
    # element_ion_pot = float(ionization_potential.loc[ionization_potential[1] == "H"][3])
#     phi = saha("H", temperature[i], partition_function, temp_arr, element_ion_pot)
    # Pe_initial = np.sqrt(phi*Pg[i])

    # #converge to a Pe value 
    # Pe_array.append(Pe_converge(Pg[i], temperature[i], Pe_initial))
    # #determine corresponding n_p value
    # n_p.append(n_H[i]*(phi/Pe_array[i]))

# #reproduce the top part of Figure 8.8
# plt.plot(VALIIIC["h"], np.log10(n_e), label = '$n_{e}$', color = 'k')
# plt.plot(VALIIIC["h"], np.log10(n_p), label = '$n_{p}$', color = 'r')
# plt.xlim(0,800)
# plt.ylim(8,14)
# plt.gca().invert_xaxis()
# plt.xlabel("h (km)")
# plt.ylabel("log n")
# plt.legend()
# plt.savefig("Figures/hw_nine_figures/fig8_8.pdf")
# plt.show()

# #ionization temperature at h = 800km
# output_array = np.stack((temperature[1:,], np.log10(n_p)[1:,]), axis=1)
# output_array=output_array[output_array[:, 0].argsort()]
# cs = CubicSpline(output_array[:, 0], output_array[:, 1])

# initial_temp = 112000
# y_diff = np.abs(cs(initial_temp) - 11)
# tol = 10**-3
# while y_diff > tol:
#     initial_temp = initial_temp + 2
#     y_diff = np.abs(cs(initial_temp) - 11)
# print(initial_temp)

"""
c) Ignoring this issue of NLTE, attempt to reproduce the bottom plot using the solar
abundances in the text file from Gray on Canvas. Use your function from Problem 17 to
calculate the LTE electron partial pressure Pe at every depth (i.e. do not use the value in
the VALIIIC table for this part) and use this to compute the ionization fractions of Fe,
Mg, Si, and H. At what heights are NLTE effects apparently important?
"""

elements_arr = ["Fe", "Mg", "Si", "H"]
n_contributions_arr = [[], [], [], []]
Pe_element_arr = [[], [], [], []]
for j in range(0, len(elements_arr)):
    for i in range(0, len(temperature)):
        element_ion_pot = float(ionization_potential.loc[ionization_potential[1] == elements_arr[j]][3])
        phi = saha(elements_arr[j], temperature[i], partition_function, temp_arr, element_ion_pot)
        Pe_initial = np.sqrt(phi*Pg[i])

        #converge to a Pe value 
        Pe_current = Pe_converge(Pg[i], temperature[i], Pe_initial)
        Pe_element_arr[j].append(Pe_current)
        #determine corresponding n_p value
        abundance = float(solar_abundance.loc[solar_abundance['element'] == elements_arr[j]]['A'])
        n_contributions_arr[j].append(abundance*n_H[i]*(phi/Pe_current))        

    plt.plot(VALIIIC["h"], np.log10(n_contributions_arr[j]/n_e), label = elements_arr[j])
plt.xlim(0,800)
plt.ylim(0,1)
plt.gca().invert_xaxis()
plt.xlabel("h (km)")
plt.ylabel("contributions to $n_{e}$")
plt.legend()
plt.savefig("Figures/hw_nine_figures/fig8_8b.pdf")
plt.show()

"""
19. Opacity and Pressure in the VALIIIC Atmosphere
a) Plot your calculated Pe from problem 18(c) versus height in the VALiiiC atmosphere.
Overplot the electron pressure implied by the perfect gas law given the electron density and
temperature listed in VALiiiC. Where are NLTE ionization effects apparently important?
"""
for j in range(0, len(elements_arr)):
    plt.plot(VALIIIC["h"], np.log10(Pe_element_arr[j]), label = elements_arr[j])
plt.plot(VALIIIC["h"], n_e*8.3145*10**7*temperature, label = 'perfect ')
plt.xlim(0,800)
plt.ylim(0,1)
plt.gca().invert_xaxis()
# plt.xlabel("h (km)")
# plt.ylabel("contributions to $n_{e}$")
plt.legend()
# plt.savefig("Figures/hw_nine_figures/fig8_8b.pdf")
plt.show()
