import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from A530_package import extinction_coefficient
from A530_package import Pe_equation, Pe_converge, saha, partition
from A530_package import opacity_total, opacity_neg_H_bf, opacity_neg_H_ff, opacity_H_ff, opacity_H_bf, opacity_electron

"""
20. Line Extinction
Calculate the monochromatic line extinction coefficient per particle for
the Na i doublet (the D line) as a function of wavelength. Include natural, quadratic Stark,
and van der Waals broadening, and assume microturbulence and thermal broadening
appropriate for the Solar photosphere. 
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
k = 1.38 * 10**-16
h = 6.62 * 10**-27

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
    first = extinction_coefficient(first_A, first_g_upper, first_g_lower, c/first_lambda, c/i, temperature, element_ion_pot, first_C4, Pe, Pg)
    second = extinction_coefficient(second_A, second_g_upper, second_g_lower, c/second_lambda, c/i, temperature, element_ion_pot, second_C4, Pe, Pg)
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

"""
21. Line and continuous opacity
Compare the opacity of line and continuous opacities for a range of wavelengths near the doublet
at the τ500 = 1 surface of the VALiiiC atmosphere. Do this by plotting the line opacity,
continuous opacities, and the sum of the continuous opacities on the same figure as a
function of wavelength. 
"""

A = solar_abundance['A'].sum()
    
#collect the values for each opacity term 
opacity_neg_H_bf_arr = []
opacity_neg_H_ff_arr = []
opacity_H_ff_arr = []
opacity_H_bf_arr = []
opacity_electron_arr = []

fe_first = (first_g_upper / 10**(partition("Na", temperature, partition_function, temp_arr)))
fe_second = (second_g_upper / 10**(partition("Na", temperature, partition_function, temp_arr)))

Na_A = float(solar_abundance.loc[solar_abundance['element'] == "Na"]['A'])
n_H = 1.166*10**17
rho = 2.727*10**-7 * n_H

continuous_opacities = []
first_line_opacities = []
second_line_opacities = []
wavelength_arr_ang = np.linspace(5800, 6000, 500) #fix 
for i in range(0, len(wavelength_arr_ang)):
    neg = False
    neg_H_bf_value = opacity_neg_H_bf(Pe, wavelength_arr_ang[i], temperature)
    if neg_H_bf_value >= 0.0 and neg == False:
        opacity_neg_H_bf_arr.append(opacity_neg_H_bf(Pe, wavelength_arr_ang[i], temperature))
    else: 
        opacity_neg_H_bf_arr.append(0.0)
        neg = True

    opacity_neg_H_ff_arr.append(opacity_neg_H_ff(Pe, wavelength_arr_ang[i], temperature))
    opacity_H_ff_arr.append(opacity_H_ff(wavelength_arr_ang[i], temperature))
    opacity_H_bf_arr.append(opacity_H_bf(wavelength_arr_ang[i], temperature))
    opacity_electron_arr.append(opacity_electron(Pe, Pg, A))

    #compute total continuous opacity from hydrogen given electron pressue, temperature, and wavelength
    opacity_total_arr = opacity_total(np.array(opacity_H_bf_arr[i]), np.array(opacity_H_ff_arr[i]), np.array(opacity_neg_H_bf_arr[i]), np.array(opacity_neg_H_ff_arr[i]), temperature, wavelength_arr_ang[i], Pe, np.array(opacity_electron_arr[i]))

    continuous_opacities.append(opacity_total_arr/(2.2701*10**(-24)))

    freq_cgs = c / (wavelength_arr_ang[i] / (10**8))
    first_line_opacities.append(((extinction_coef_first[i]*Na_A*n_H*fe_first*(phi/Pe))/rho)*(1-np.exp(-(h*freq_cgs)/(k*temperature))))
    second_line_opacities.append(((extinction_coef_second[i]*Na_A*n_H*fe_second*(phi/Pe))/rho)*(1-np.exp(-(h*freq_cgs)/(k*temperature))))

plt.plot(wavelength_arr_ang, first_line_opacities, linewidth = 2, label = 'λ = 5889.95', color = 'r')
plt.plot(wavelength_arr_ang, second_line_opacities, linewidth = 2, label = 'λ = 5895.924', color = 'g')
plt.plot(wavelength_arr_ang, continuous_opacities, linewidth = 2, label = 'continuous', color = 'k')
plt.xlabel("wavelength [Å]", fontsize = 12)
plt.ylabel("opacity", fontsize = 12)
plt.xticks(fontsize = 12) 
plt.yticks(fontsize = 12) 
plt.legend(fontsize = 12)
plt.legend()
plt.yscale('log')
plt.savefig("Figures/hw_ten_figures/opacity.pdf", bbox_inches='tight')
plt.show()