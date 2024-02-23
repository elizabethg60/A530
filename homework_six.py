import numpy as np
import pandas as pd
from A530_package import partition, saha

"""
11. Partition Functions
Write code to calculate the partition function for an arbitrary species at an arbitrary
temperature. 
"""

#read in data
partition_function = pd.read_csv('data/RepairedPartitionFunctions.txt', sep= ' ', header=None)

#create temp array
theta_arr = np.arange(0.2, 2.2, step = 0.2)
temp_arr = 5040/theta_arr

#calculate the partition function for an arbitrary species at an arbitrary temperature
element = "H"
temperature = 4000
# #returns log of the partition function
# print("log partition function: {}".format(partition(element, temperature, partition_function, temp_arr)))
# print("partition function: {}".format(10**(partition(element, temperature, partition_function, temp_arr))))

"""
13. Saha's Φ(T)
Write a function Φ(T) that returns the temperature dependent part of the RHS of the
Saha equation for an arbitrary element.
"""

#read in data 
ionization_potential = pd.read_fwf('data/ioniz.txt', header=None)

#determine (first) ionization potential for given element (no second / third ionization potential at this time)
if element == "H-":
    element_ion_pot = 0.755
else: 
    element_ion_pot = float(ionization_potential.loc[ionization_potential[1] == element][3])

#return temperature dependent part of RHS of Saha equation 
print("Saha temperature dependent: {}".format(saha(element, temperature, partition_function, temp_arr, element_ion_pot)))