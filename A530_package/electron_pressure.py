import numpy as np
import pandas as pd
from A530_package import saha

path = '/Users/efg5335/Desktop/Courses/A530/data/'

def Pe_equation(Pg, temperature, Pe):
# returns the RHS of equation 9.9 in Gray given a Pg, temp, and Pe

    #read in data
    partition_function = pd.read_csv(path + 'RepairedPartitionFunctions.txt', sep= ' ', header=None)
    ionization_potential = pd.read_fwf(path + 'ioniz.txt', header=None)
    solar_abundance = pd.read_csv(path + 'Solar_Data/SolarAbundance.txt', sep="	")

    #create temp array for phi function 
    theta_arr = np.arange(0.2, 2.2, step = 0.2)
    temp_arr = 5040/theta_arr

    #create element array of elements needed in Pe calculation 
    element_arr = ["H", "He", "C", "Si", "Fe", "Mg", "Ni", "Cr", "Ca", "Na", "K"]

    #collect components of RHS for each element 
    top = []
    bottom = []
    for element in element_arr: 
        #determine ionization potential for given element
        element_ion_pot = float(ionization_potential.loc[ionization_potential[1] == element][3])
        #determine solar abundance for given element
        A = float(solar_abundance.loc[solar_abundance["element"] == element]["A"])

        #compute saha phi value for given inputs 
        phi = saha(element, temperature, partition_function, temp_arr, element_ion_pot)

        top.append((A*phi/Pe)/(1+phi/Pe))
        bottom.append(A*(1+ (phi/Pe)/(1+phi/Pe)))

    return Pg*(np.sum(top)/np.sum(bottom))

def Pe_converge(Pg, temperature, Pe):
# returns the converged Pe value using iterations of equation 9.9 in Gray

    #set convergence tolerance
    tol = 10**(-16)
    #initial difference between LHS and RHS of equation 9.9
    diff = np.abs(Pe - Pe_equation(Pg, temperature, Pe))
    #iterate till difference meets the convergence tolerance
    while diff > tol:
        #update estimated Pe value as mean of previous LHS and RHS values 
        Pe = np.mean([Pe, Pe_equation(Pg, temperature, Pe)])
        diff = np.abs(Pe - Pe_equation(Pg, temperature, Pe))
    
    return Pe