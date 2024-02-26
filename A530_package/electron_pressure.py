import numpy as np
import pandas as pd
from A530_package import partition, saha

path = '/Users/efg5335/Desktop/Courses/A530/data/'

def Pe_equation(Pg, temperature, Pe):
    #read in data
    partition_function = pd.read_csv(path + 'RepairedPartitionFunctions.txt', sep= ' ', header=None)
    ionization_potential = pd.read_fwf(path + 'ioniz.txt', header=None)
    solarabundance = pd.read_csv(path + 'Solar_Data/SolarAbundance.txt', sep="	")

    #create temp array
    theta_arr = np.arange(0.2, 2.2, step = 0.2)
    temp_arr = 5040/theta_arr

    element_arr = ["H", "He", "C", "Si", "Fe", "Mg", "Ni", "Cr", "Ca", "Na", "K"]

    top = []
    bottom = []
    for element in element_arr: 
        #determine (first) ionization potential for given element (no second / third ionization potential at this time)
        if element == "H-":
            element_ion_pot = 0.755
        else: 
            element_ion_pot = float(ionization_potential.loc[ionization_potential[1] == element][3])

        phi = saha(element, temperature, partition_function, temp_arr, element_ion_pot)
        A = float(solarabundance.loc[solarabundance["element"] == element]["A"])

        top.append((A*phi/Pe)/(1+phi/Pe))
        bottom.append(A*(1+ (phi/Pe)/(1+phi/Pe)))

    return Pg*(np.sum(top)/np.sum(bottom))