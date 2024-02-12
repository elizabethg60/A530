import pandas as pd

"""
11. Partition Functions
Write code to calculate the partition function for an arbitrary species at an arbitrary
temperature. We will rarely deal with double ionization, so your code does not have to
successfully handle every case of a doubly ionized element. Hints:
• When you solve the Saha equation for the ionization levels of atoms, especially
hydrogen, you will need to have values for partition functions of things like Hii, and
H-, so be sure your partition function returns useful numbers for these.
• H- has an ionization potential of 0.755 eV and you will recall that it has exactly
one bound state — what do you suppose its electron configuration is?
• What does this imply about the occupation number g0 and the partition function
UH- (T )?
"""

partition_function = pd.read_csv('data/RepairedPartitionFunctions.txt', sep= ' ', header=None)

#to do: (1) confirm with Gray's appendix and understand parition table or is it statistical weights 
#       (2) understand problem (check Gray, Rutten is pg 30 only)
#       (3) work on problem THEN consider hints 

"""
13. Saha's Φ(T)
Write a function Φ(T) that returns the temperature dependent part of the RHS of the
Saha equation for an arbitrary element.
"""

ionization_potential = pd.read_fwf('data/ioniz.txt', header=None)
#Atomic number, Element, Average atomic weight, 1st, 2nd, and 3rd ionization energies in eV.  
#to get specific element index by row number - 1 using .iloc
#ionization_potential.loc[ionization_potential[1] == "element"]

#to do: (1) solve saha equation for temperature (2) figure out LHS (3) write function and finish answering question