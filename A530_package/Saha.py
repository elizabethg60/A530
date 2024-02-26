import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

def partition(element, temperature, data, temp_arr):
# returns partition function given an arbitrary element and temperature 

    #HII case
    if element == "HII":
        return 0.0
    #H- case
    if element == "H-":
        return -(5040/temperature)*(0.755)

    #if element not found, return log partition function of 0.0
    if len(data.loc[data[0] == element]) == 0.0:
        return 0.0

    #collect true parition functions from table
    true_partition = np.array(data.loc[data[0] == element])[0][1:-1]

    partition_arr = []
    start = 0
    for i in range(0, len(true_partition)):
        if true_partition[i] != "-":
            partition_arr.append(float(true_partition[i]))
        if true_partition[i] == "-":
            start = i + 1

    #determine new array of temperatures for known partitions 
    new_temp_arr = temp_arr[start:,]

    #fit a cubic spline to data to get partition function at a given temp
    cs = CubicSpline(new_temp_arr[::-1], partition_arr[::-1])
    #plot for confirmation 
    # plt.scatter(new_temp_arr, partition_arr)
    # x_range = np.arange(min(new_temp_arr), max(new_temp_arr), 0.1)
    # plt.plot(x_range, cs(x_range), label='Cubic Spline')
    # plt.xlabel("temperature (K)")
    # plt.ylabel("log partition function")
    # plt.show()

    return cs(temperature)

def saha(element, temperature, data, temp_arr, ion_potential):
# returns temperature dependent part of RHS of Saha equation given an arbitrary element and temperature

    #H- case
    if element == "H-":
        U = 10**partition(element, temperature, data, temp_arr)
        U_plus = 10**partition("H", temperature, data, temp_arr)

    #every other case
    else: 
        U = 10**partition(element, temperature, data, temp_arr)
        U_plus = 10**partition(element+"+", temperature, data, temp_arr) 

    return 0.6665*(U_plus/U)*(temperature**(5/2))*10**(-(5040/temperature)*ion_potential)