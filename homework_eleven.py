import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from A530_package import extinction_coefficient, tau_eta
from A530_package import Pe_equation, Pe_converge, saha, partition
from A530_package import quadratic_source, linear_source, eddington_flux
from A530_package import opacity_total, opacity_neg_H_bf, opacity_neg_H_ff, opacity_H_ff, opacity_H_bf, opacity_electron 

"""
22. The Sodium Doublet in LTE
Calculate the emergent flux in the VALiiiC atmosphere for the Na i doublet. Use the
VALiiiC value for ne (and so Pe from the perfect gas law). Use your work from Problem
7, your continuous opacity, and your sodium opacity from previous problems to calculate
τν (τ500, v) and F +v (0). Show a well-chosen plot for τν vs. τ500 and v and plot the emergent
flux vs. λ. Does your answer look like the Solar D doublet? (Rutten shows what the
doublet looks like).
The VALiiiC atmosphere only has 50 or so grid points. Youll need to build a finer grid
for this problem via interpolation. Where these points are needed will depend on the
frequency youre interested in. I added 1000 points between 10-8 and 10-2.
"""

#read in data
path = '/Users/efg5335/Desktop/Courses/A530/data/'
solar_abundance = pd.read_csv(path + 'Solar_Data/SolarAbundance.txt', sep="	")
Na_A = float(solar_abundance.loc[solar_abundance['element'] == "Na"]['A'])
A = solar_abundance['A'].sum()

partition_function = pd.read_csv(path + 'RepairedPartitionFunctions.txt', sep= ' ', header=None)
ionization_potential = pd.read_fwf(path + 'ioniz.txt', header=None)
element_ion_pot = float(ionization_potential.loc[ionization_potential[1] == "Na"][3])

theta_arr = np.arange(0.2, 2.2, step = 0.2)
temp_arr = 5040/theta_arr

VALIIIC = pd.read_csv(path + 'Solar_Data/VALIIIC_sci_e.txt', sep= ' ', header=None)
VALIIIC.set_axis(["h", "m", "tau_500", "T", "V", "n_H", "n_e", "Ptotal", "Pgas/Ptotal", "rho"], axis=1, inplace=True)
tau_500 = VALIIIC["tau_500"]

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

#hardwire values for a_n coefficients
a0 = a1 = 1
a2 = 2

wavelength_range = np.linspace(5.8*10**-5, 6*10**-5, 500)

tau_vs_tau500 = []
# for ind in range(0, len(tau_500)):#[48]:
#     #hardwire temperature and gas pressure 
#     temperature = VALIIIC["T"][ind]
#     Pg = VALIIIC["Ptotal"][ind] * VALIIIC["Pgas/Ptotal"][ind]

#     #derive converged electron pressure 
#     phi = saha("Na", temperature, partition_function, temp_arr, element_ion_pot)
#     Pe_initial = np.sqrt(phi*Pg)
#     Pe = Pe_converge(Pg, temperature, Pe_initial)

#     #iterate through a range of wavelength and determine extinction coefficient 
#     extinction_coef_first = []
#     extinction_coef_second = []
#     continuous_opacities = []
#     for i in [wavelength_range[0]]:
#         first = extinction_coefficient(first_A, first_g_upper, first_g_lower, c/first_lambda, c/i, temperature, element_ion_pot, first_C4, Pe, Pg)
#         second = extinction_coefficient(second_A, second_g_upper, second_g_lower, c/second_lambda, c/i, temperature, element_ion_pot, second_C4, Pe, Pg)
#         extinction_coef_first.append(first)
#         extinction_coef_second.append(second)

#     #collect the values for each opacity term 
#     fe_ground = (second_g_upper / 10**(partition("Na", temperature, partition_function, temp_arr)))
#     n_H = VALIIIC["n_H"][ind]
#     rho = VALIIIC["rho"][ind]

#     first_line_opacities = []
#     second_line_opacities = []
#     continuous_opacities = []

#     opacity_neg_H_bf_arr = []
#     opacity_neg_H_ff_arr = []
#     opacity_H_ff_arr = []
#     opacity_H_bf_arr = []
#     opacity_electron_arr = []
#     for i in range(0, len([wavelength_range[0]])):
#         lambda_ang = wavelength_range[i]*10**8

#         neg = False
#         neg_H_bf_value = opacity_neg_H_bf(Pe, lambda_ang, temperature)
#         if neg_H_bf_value >= 0.0 and neg == False:
#             opacity_neg_H_bf_arr.append(opacity_neg_H_bf(Pe, lambda_ang, temperature))
#         else: 
#             opacity_neg_H_bf_arr.append(0.0)
#             neg = True

#         opacity_neg_H_ff_arr.append(opacity_neg_H_ff(Pe, lambda_ang, temperature))
#         opacity_H_ff_arr.append(opacity_H_ff(lambda_ang, temperature))
#         opacity_H_bf_arr.append(opacity_H_bf(lambda_ang, temperature))
#         opacity_electron_arr.append(opacity_electron(Pe, Pg, A))

#         #compute total continuous opacity from hydrogen given electron pressue, temperature, and wavelength
#         opacity_total_arr = opacity_total(np.array(opacity_H_bf_arr[i]), np.array(opacity_H_ff_arr[i]), np.array(opacity_neg_H_bf_arr[i]), np.array(opacity_neg_H_ff_arr[i]), temperature, lambda_ang, Pe, np.array(opacity_electron_arr[i]))

#         continuous_opacities.append(opacity_total_arr/(2.2701*10**(-24)))

#         freq_cgs = c / wavelength_range[i]
#         first_line_opacities.append(((extinction_coef_first[i]*Na_A*n_H*fe_ground*(1/(1+(phi/Pe))))/rho)*(1-np.exp(-(h*freq_cgs)/(k*temperature))))
#         second_line_opacities.append(((extinction_coef_second[i]*Na_A*n_H*fe_ground*(1/(1+(phi/Pe))))/rho)*(1-np.exp(-(h*freq_cgs)/(k*temperature))))

#     eta = (np.array(first_line_opacities) + np.array(second_line_opacities)) / continuous_opacities
#     tau_array = tau_eta(eta)
#     tau_vs_tau500.append(tau_array)
#     # emergent_flux = eddington_flux(quadratic_source, tau_array, a0, a1, a2) * 4 * np.pi

#     # plt.plot(wavelength_range*10**8, tau_array, linewidth = 2, color = 'k')
#     # plt.xlabel("wavelength [Å]", fontsize = 12)
#     # plt.ylabel("tau", fontsize = 12)
#     # plt.xticks(fontsize = 12) 
#     # plt.yticks(fontsize = 12) 
#     # plt.yscale('log')
#     # plt.savefig("Figures/hw_eleven_figures/tau_v_lambda.pdf", bbox_inches='tight')
#     # plt.show()

# plt.plot(tau_500, tau_vs_tau500, linewidth = 2, color = 'k')
# plt.xlabel("τ500", fontsize = 12)
# plt.ylabel("tau", fontsize = 12)
# plt.xticks(fontsize = 12) 
# plt.yticks(fontsize = 12) 
# plt.title('wavelength: {} [Å]'.format(wavelength_range[0]*10**8))
# plt.savefig("Figures/hw_eleven_figures/tau_v_tau500.pdf", bbox_inches='tight')
# plt.show()

emergent_flux_arr = []
for i in range(0, len(wavelength_range)):
    extinction_coef_first = []
    extinction_coef_second = []
    continuous_opacities = []

    first_line_opacities = []
    second_line_opacities = []
    continuous_opacities = []

    opacity_neg_H_bf_arr = []
    opacity_neg_H_ff_arr = []
    opacity_H_ff_arr = []
    opacity_H_bf_arr = []
    opacity_electron_arr = []

    for ind in range(0, len(tau_500)):

        #hardwire temperature and gas pressure 
        temperature = VALIIIC["T"][ind]
        Pg = VALIIIC["Ptotal"][ind] * VALIIIC["Pgas/Ptotal"][ind]
        n_e = VALIIIC["n_e"][ind]

        #derive converged electron pressure 
        phi = saha("Na", temperature, partition_function, temp_arr, element_ion_pot)
        Pe = n_e*1.38*10**(-16)*temperature

        #iterate through a range of wavelength and determine extinction coefficient 
        first = extinction_coefficient(first_A, first_g_upper, first_g_lower, c/first_lambda, c/wavelength_range[i], temperature, element_ion_pot, first_C4, Pe, Pg)
        second = extinction_coefficient(second_A, second_g_upper, second_g_lower, c/second_lambda, c/wavelength_range[i], temperature, element_ion_pot, second_C4, Pe, Pg)
        extinction_coef_first.append(first)
        extinction_coef_second.append(second)    

        #collect the values for each opacity term 
        fe_ground = (second_g_upper / 10**(partition("Na", temperature, partition_function, temp_arr)))
        n_H = VALIIIC["n_H"][ind]
        rho = VALIIIC["rho"][ind]

        lambda_ang = wavelength_range[i]*10**8

        neg = False
        neg_H_bf_value = opacity_neg_H_bf(Pe, lambda_ang, temperature)
        if neg_H_bf_value >= 0.0 and neg == False:
            opacity_neg_H_bf_arr.append(opacity_neg_H_bf(Pe, lambda_ang, temperature))
        else: 
            opacity_neg_H_bf_arr.append(0.0)
            neg = True

        opacity_neg_H_ff_arr.append(opacity_neg_H_ff(Pe, lambda_ang, temperature))
        opacity_H_ff_arr.append(opacity_H_ff(lambda_ang, temperature))
        opacity_H_bf_arr.append(opacity_H_bf(lambda_ang, temperature))
        opacity_electron_arr.append(opacity_electron(Pe, Pg, A))

        #compute total continuous opacity from hydrogen given electron pressue, temperature, and wavelength
        opacity_total_arr = opacity_total(np.array(opacity_H_bf_arr[ind]), np.array(opacity_H_ff_arr[ind]), np.array(opacity_neg_H_bf_arr[ind]), np.array(opacity_neg_H_ff_arr[ind]), temperature, lambda_ang, Pe, np.array(opacity_electron_arr[ind]))

        continuous_opacities.append(opacity_total_arr/(2.2701*10**(-24)))

        freq_cgs = c / wavelength_range[i]
        first_line_opacities.append(((extinction_coef_first[ind]*Na_A*n_H*fe_ground*(1/(1+(phi/Pe))))/rho)*(1-np.exp(-(h*freq_cgs)/(k*temperature))))
        second_line_opacities.append(((extinction_coef_second[ind]*Na_A*n_H*fe_ground*(1/(1+(phi/Pe))))/rho)*(1-np.exp(-(h*freq_cgs)/(k*temperature))))

    eta = (np.array(first_line_opacities) + np.array(second_line_opacities)) / continuous_opacities
    tau_array = tau_eta(eta)
    emergent_flux_arr.append(eddington_flux(quadratic_source, tau_array[1:,], a0, a1, a2) * 4 * np.pi)

plt.plot(wavelength_range*10**8, emergent_flux_arr, linewidth = 2, color = 'k')
plt.xlabel("wavelength [Å]", fontsize = 12)
plt.ylabel("emergent flux", fontsize = 12)
plt.xticks(fontsize = 12) 
plt.yticks(fontsize = 12) 
plt.savefig("Figures/hw_eleven_figures/emer_flux_2.pdf", bbox_inches='tight')
plt.show()