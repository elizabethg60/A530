import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from A530_package import Pe_converge, saha
from A530_package import opacity_total, opacity_neg_H_bf, opacity_neg_H_ff, opacity_H_ff, opacity_H_bf, opacity_electron

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

Pe_array = []
n_p = []
for i in range(0, len(temperature)):
    #initial guess for Pe 
    element_ion_pot = float(ionization_potential.loc[ionization_potential[1] == "H"][3])
    phi = saha("H", temperature[i], partition_function, temp_arr, element_ion_pot)
    Pe_initial = np.sqrt(phi*Pg[i])

    #converge to a Pe value 
    Pe_array.append(Pe_converge(Pg[i], temperature[i], Pe_initial))
    #determine corresponding n_p value
    n_p.append(n_H[i]*(phi/Pe_array[i]))

#reproduce the top part of Figure 8.8
plt.plot(VALIIIC["h"], np.log10(n_e), label = '$n_{e}$', color = 'k')
plt.plot(VALIIIC["h"], np.log10(n_p), label = '$n_{p}$', color = 'r')
plt.xlim(0,800)
plt.ylim(8,14)
plt.gca().invert_xaxis()
plt.xlabel("h (km)")
plt.ylabel("log n")
plt.legend()
plt.savefig("Figures/hw_nine_figures/fig8_8.pdf")
plt.show()

#ionization temperature at h = 800km
output_array = np.stack((temperature[1:,], np.log10(n_p)[1:,]), axis=1)
output_array=output_array[output_array[:, 0].argsort()]
cs = CubicSpline(output_array[:, 0], output_array[:, 1])

initial_temp = 112000
y_diff = np.abs(cs(initial_temp) - 11)
tol = 10**-3
while y_diff > tol:
    initial_temp = initial_temp + 2
    y_diff = np.abs(cs(initial_temp) - 11)
print(initial_temp)

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
        n_contributions_arr[j].append(abundance*n_H[i]*((phi/Pe_current)/(1+(phi/Pe_current))))        

    plt.plot(VALIIIC["h"], np.log10(n_contributions_arr[j]/n_e), label = elements_arr[j])
plt.xlim(0,800)
plt.ylim(-1,1)
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
Pe_sum_arr = [sum(i) for i in zip(*Pe_element_arr)]
plt.plot(VALIIIC["h"], np.log10(Pe_sum_arr), label = 'computed')
plt.plot(VALIIIC["h"], np.log10(n_e*1.38*10**(-16)*temperature), label = 'perfect gas law')
plt.xlim(0,800)
plt.ylim(-5,5)
plt.gca().invert_xaxis()
plt.xlabel("h (km)")
plt.ylabel("$P_{e}$")
plt.legend()
plt.savefig("Figures/hw_nine_figures/fig19a.pdf")
plt.show()

"""
b) Combine your work from problems 15 and 17 into a function that calculates the total
continuous opacity (in its proper units, i.e. cm2/g) including electron scattering, given a
wavelength (or frequency), vector of abundances, gas pressure, and electron pressure.

c) Check your new opacity function against the numbers in the following table for con-
ditions on row 48 of the VALiiiC atmosphere and for one of the plots from Problem 15.
In these opacities, I have included all ionization and stimulated emission factors, so the
total (“continuum”) opacity is the sum of the numbers below it. If you feel that I must
have miscalculated a number, please let me know which one and what factors go into that
number. You can also confirm that your electron pressure routine works from the gas
pressure in the second column
"""

A = 0
A_weighted = (solar_abundance['A'] * solar_abundance['weight']).sum()

for j in range(0, len(elements_arr)):
    A = A + float(solar_abundance.loc[solar_abundance['element'] == elements_arr[j]]['A'])
    
#hardwire temperature and electron pressure according to table in hw 9
temp = 6420
Pe = 57.0
Pg_single = 1.13*10**5

#collect the values for each opacity term 
opacity_neg_H_bf_arr = []
opacity_neg_H_ff_arr = []
opacity_H_ff_arr = []
opacity_H_bf_arr = []
opacity_electron_arr = []

wavelength = 5000
neg = False
neg_H_bf_value = opacity_neg_H_bf(Pe, wavelength, temp)
if neg_H_bf_value >= 0.0 and neg == False:
    opacity_neg_H_bf_arr.append(opacity_neg_H_bf(Pe, wavelength, temp))
else: 
    opacity_neg_H_bf_arr.append(0.0)
    neg = True

opacity_neg_H_ff_arr.append(opacity_neg_H_ff(Pe, wavelength, temp))
opacity_H_ff_arr.append(opacity_H_ff(wavelength, temp))
opacity_H_bf_arr.append(opacity_H_bf(wavelength, temp))
opacity_electron_arr.append(opacity_electron(Pe, Pg_single, A))

#compute total continuous opacity from hydrogen given electron pressue, temperature, and wavelength
opacity_total_arr = opacity_total(np.array(opacity_H_bf_arr), np.array(opacity_H_ff_arr), np.array(opacity_neg_H_bf_arr), np.array(opacity_neg_H_ff_arr), temp, wavelength, Pe, np.array(opacity_electron_arr))

print(opacity_total_arr/(2.2701*10**(-24)))

"""
d) Using your Pe values from Problem 18(c), calculate κ500 as a function of τ500 with your
now-optimized opacity calculator.

e) Check your answer by inferring κ500 from the run of Ptotal and τ500 in the VALiiiC atmo-
sphere and the equation of hydrostatic equilibrium. Show how your results compare; why
are they not in perfect agreement? (Hints: In the VALiiiC atmosphere, log g = 4.4377.
Determine dτ500 as the step in τ500 between rows in the VALiiiC atmosphere tabulation,
and determine the state variables at those points by averaging their neighboring values.
That is, calculate dP , dτ500, T , and Ptotal in 51 bins where the value for T in the first bin
is the average of T in the first two rows tabulation.)
"""

κ500 = []
for i in range(0, len(temperature)):
    temp = temperature[i]
    Pe = Pe_sum_arr[i]
    Pg_current = Pg[i]

    #collect the values for each opacity term 
    opacity_neg_H_bf_arr = []
    opacity_neg_H_ff_arr = []
    opacity_H_ff_arr = []
    opacity_H_bf_arr = []
    opacity_electron_arr = []

    wavelength = 5000
    neg = False
    neg_H_bf_value = opacity_neg_H_bf(Pe, wavelength, temp)
    if neg_H_bf_value >= 0.0 and neg == False:
        opacity_neg_H_bf_arr.append(opacity_neg_H_bf(Pe, wavelength, temp))
    else: 
        opacity_neg_H_bf_arr.append(0.0)
        neg = True

    opacity_neg_H_ff_arr.append(opacity_neg_H_ff(Pe, wavelength, temp))
    opacity_H_ff_arr.append(opacity_H_ff(wavelength, temp))
    opacity_H_bf_arr.append(opacity_H_bf(wavelength, temp))
    opacity_electron_arr.append(opacity_electron(Pe, Pg_current, A))

    #compute total continuous opacity from hydrogen given electron pressue, temperature, and wavelength
    opacity_total_arr = opacity_total(np.array(opacity_H_bf_arr), np.array(opacity_H_ff_arr), np.array(opacity_neg_H_bf_arr), np.array(opacity_neg_H_ff_arr), temp, wavelength, Pe, np.array(opacity_electron_arr))
    κ500.append(opacity_total_arr/A_weighted)

plt.plot(VALIIIC["tau_500"], np.log10(κ500))

dPdtau = np.diff(VALIIIC["Ptotal"])/np.diff(VALIIIC["tau_500"])
g = 10**4.4377

plt.plot(VALIIIC["tau_500"][1:,], np.log10(g/dPdtau))
plt.xlabel("τ500")
plt.ylabel("κ500")
plt.savefig("Figures/hw_nine_figures/fig19d.pdf")
plt.show()