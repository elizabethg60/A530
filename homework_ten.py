import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from A530_package import extinction_coefficient
from A530_package import Pe_equation, Pe_converge, saha

"""
20. Line Extinction
Calculate the monochromatic line extinction coefficient per particle for
the Na i doublet (the D line) as a function of wavelength. Include natural, quadratic Stark,
and van der Waals broadening, and assume microturbulence and thermal broadening
appropriate for the Solar photosphere. Be sure to document your work here—how did 
you calculate the various terms that contribute to the extinction coefficient? 
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

#required constants (cgs)
first_lambda = 5889.95 / (10**8)
first_g_lower = 2
first_g_upper = 4
first_A = 6.16 * 10**-1 * 10**8
first_C4 = -15.17

second_lambda = 5895.924 / (10**8)
second_g_lower = 2
second_g_upper = 2
second_A = 6.14 * 10**-1 * 10**8
second_C4 = -15.33

c = 3 * 10**10

#hardwire temperature and gas pressure 
temperature = 6420
Pg = 1.172*10**5 * 0.9702

#derive converged electron pressure 
element_ion_pot = float(ionization_potential.loc[ionization_potential[1] == "Na"][3])
phi = saha("Na", temperature, partition_function, temp_arr, element_ion_pot)
Pe_initial = np.sqrt(phi*Pg)
Pe = Pe_converge(Pg, temperature, Pe_initial)

#iterate through a range of wavelength and determine extinction coefficient 
extinction_coef_first = []
extinction_coef_second = []
extinction_coef_sum = []
wavelength_range = np.linspace(5.8*10**-5, 6*10**-5, 500)
for i in wavelength_range:
    first = extinction_coefficient(first_A, first_g_lower, first_g_upper, c/first_lambda, c/i, temperature, element_ion_pot, first_C4, Pe, Pg)
    second = extinction_coefficient(second_A, second_g_lower, second_g_upper, c/second_lambda, c/i, temperature, element_ion_pot, second_C4, Pe, Pg)
    extinction_coef_first.append(first)
    extinction_coef_second.append(second)
    extinction_coef_sum.append(first + second)

plt.plot(wavelength_range*10**8, extinction_coef_first, linewidth = 2, label = 'λ = 5889.95', color = 'r')
plt.plot(wavelength_range*10**8, extinction_coef_second, linewidth = 2, label = 'λ = 5895.924', color = 'g')
plt.plot(wavelength_range*10**8, extinction_coef_sum, linewidth = 2, label = 'Sum', color = 'k')
plt.xlabel("wavelength [Å]", fontsize = 12)
plt.ylabel("log extinction coefficient", fontsize = 12)
plt.xticks(fontsize = 12) 
plt.yticks(fontsize = 12) 
plt.legend(fontsize = 12)
plt.legend()
plt.yscale('log')
plt.savefig("Figures/hw_ten_figures/extinction.pdf", bbox_inches='tight')
plt.show()