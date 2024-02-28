import numpy as np
import pandas as pd
from A530_package import Pe_equation, Pe_converge, saha

"""
Write a program to implement Equation (9.8) from D. F. Gray (3rd edition) for parts
of the Sun. Consider the all elements for which you have partition functions, ionizations,
and abundances.
Use Table 9.2 in Gray to check your code.
"""

#hardwire temperatue and Pg
temperature = 6989
Pg = 10**(5.14)

#initial guess for Pe 
path = '/Users/efg5335/Desktop/Courses/A530/data/'
partition_function = pd.read_csv(path + 'RepairedPartitionFunctions.txt', sep= ' ', header=None)
ionization_potential = pd.read_fwf(path + 'ioniz.txt', header=None)
element_ion_pot = float(ionization_potential.loc[ionization_potential[1] == "H"][3])
theta_arr = np.arange(0.2, 2.2, step = 0.2)
temp_arr = 5040/theta_arr
phi = saha("H", temperature, partition_function, temp_arr, element_ion_pot)
Pe_initial = np.sqrt(phi*Pg)

Pe = Pe_converge(Pg, temperature, Pe_initial)

print("converged Pe value: {}".format(Pe))
print("converged logPe value: {}".format(np.log10(Pe)))

"""
What are the sums of abundances ΣAj and abundance-weighted masses ΣAj μj for the
lightest 30 elements?
"""

solar_abundance = pd.read_csv(path + 'Solar_Data/SolarAbundance.txt', sep="	")
lightest = solar_abundance.loc[solar_abundance['atomic'] <= 30]

print("sums of abundances of lightest 30 elements: {}".format(np.sum(lightest["A"])))
print("sums of abundance-weighted mass of lightest 30 elements: {}".format((lightest['A'] * lightest['weight']).sum()))

"""
What is ΣAj for the “metals” only (Z > 2)? What is ΣAj μj for the metals?
"""

metals = solar_abundance.loc[solar_abundance['atomic'] > 2]

print("sums of abundances of metals: {}".format(np.sum(metals["A"])))
print("sums of abundance-weighted mass of metals: {}".format((metals['A'] * metals['weight']).sum()))