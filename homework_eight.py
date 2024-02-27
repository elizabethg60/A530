import numpy as np
import pandas as pd
from A530_package import Pe_equation, Pe_converge

"""
Write a program to implement Equation (9.8) from D. F. Gray (3rd edition) for parts
of the Sun. Consider the all elements for which you have partition functions, ionizations,
and abundances.
Use Table 9.2 in Gray to check your code.
"""

#hardwire temperatue and Pg
temperature = 4000
Pg = 10**(5.27)

#initial guess for Pe found using equation 9.11 in Gray 
path = '/Users/efg5335/Desktop/Courses/A530/data/'
solar_abundance = pd.read_csv(path + 'Solar_Data/SolarAbundance.txt', sep="	")
element_arr = ["C", "Si", "Fe", "Mg", "Ni", "Cr", "Ca", "Na", "K"]
A = [float(solar_abundance.loc[solar_abundance["element"] == element]["A"]) for element in element_arr]
Pe_initial = (Pg*np.sum(A))/1.085

Pe = Pe_converge(Pg, temperature, Pe_initial)

print("converged Pe value: {}".format(Pe))

print(np.log10(Pe))
