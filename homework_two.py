from A530_package import Planck, box_integrator, integrate_Planck
import matplotlib.pyplot as plt
from astropy import constants
import numpy as np

h = constants.h.cgs.value # cm^2 g s^-1
c = constants.c.cgs.value # cm/s
k = constants.k_B.cgs.value # cm^2 g s^-2 K^-1

"""
(a) Numerically integrate Bν (T = 7500 K) over the interval 0 <  ̃ν < ∞. Your result will have the units of specific intensity 
times wavenumber. For this exercise, use a box or trapezoid rule integrator. Your box or trapezoid function should accept an 
array of x and an array of y points, and should not assume the points are evenly spaced (but you can assume they’re ordered in
x). Then, wrap this integrator in another function that lets you specify the the minimum and maximum x-axis value and density 
of points in the integration of a given function.
"""

#set temperature (Kelvin)
temp = 7500
#set wavenumber range: 0 to "inf"
wavenumber_microns = np.linspace(0, 1e5, int(1e7))
#wavelength range 
wavelengths = 1 / wavenumber_microns
intensity_function = Planck(c / (wavelengths/1e4), temp)
plt.plot(wavenumber_microns, intensity_function * 1e4, label = "{}K".format(temp), linewidth = 5)
plt.xlabel("wavenumber [$μm^{−1}$]", fontsize = 12)
plt.ylabel("intensity [$10^{-4}$ erg $cm^{-2}$ $s^{-1}$ $st^{-1}$ $Hz^{-1}$]", fontsize = 12)
plt.xticks(fontsize = 12) 
plt.yticks(fontsize = 12) 
plt.legend(fontsize = 12)
plt.xscale("log")
plt.savefig("Figures/hw_two_figures/sample_a.png", bbox_inches='tight')
plt.show()

#area = 0.00019049615513313137 (units of specific intensity times wavenumber)
print("area under curve: {} (specific intensity times wavenumber)".format(box_integrator(wavenumber_microns, intensity_function))) 
print("area under curve: {} (specific intensity times wavenumber)".format(integrate_Planck(temp, 0, 1e5, int(1e7)))) 

"""
b) Determine the precision of your calculation through comparison with the analytic result for the area under the Planck function. 
Precision is measured as |(truth−calculation)|/truth =|1 − calculation/truth|, and is most usefully plotted logarithmically.
Your integrator has three parameters you can vary: which one(s) is(are) limiting your precision? Can you get to machine precision? 
Plot your precision vs. the 3 parameters with the others held constant. Interpret these.
"""

#using A502 homework one notes + Christian's help for unit conversion into specific intensity times wavenumber
truth = (2*h*(k**4)*(np.pi**4)*(temp**4))/(15*(c**2)*(h**4))/c/10**4 

calculation = integrate_Planck(temp, 0, 1e5, int(1e7))
precision = np.abs(truth - calculation)/truth
print("truth: {} vs calculation: {} gives precision of {}".format(truth, calculation, precision))
#truth: 0.000190496155149685 vs calculation: 0.00019049615513313137 gives precision of 8.689742772790219e-11

x_min_arr = np.linspace(0, 0.1, 20)
x_max_arr = np.linspace(1e5, 1e7, 20)
density_arr = np.linspace(1e7, 1e8, 20)

x_min_result = []
x_max_result = []
density_result = []
for i in x_min_arr:
    calculation = integrate_Planck(temp, i, 1e5, int(1e7))
    precision = np.abs(truth - calculation)/truth
    x_min_result.append(precision)
for i in x_max_arr:
    calculation = integrate_Planck(temp, 0, i, int(1e7))
    precision = np.abs(truth - calculation)/truth
    x_max_result.append(precision)
for i in density_arr:
    calculation = integrate_Planck(temp, 0, 1e5, int(i))
    precision = np.abs(truth - calculation)/truth
    density_result.append(precision)

#precision worsens as min wavenumber increases
plt.plot(x_min_arr, x_min_result, linewidth = 1)
plt.scatter(x_min_arr, x_min_result, s = 12)
plt.xlabel("min wavenumber [$μm^{−1}$]", fontsize = 12)
plt.ylabel("log precision", fontsize = 12)
plt.xticks(fontsize = 12) 
plt.yticks(fontsize = 12) 
plt.yscale("log")
plt.savefig("Figures/hw_two_figures/sample_bi.png", bbox_inches='tight')
plt.show()

#precision worsens as max wavenumber increases (given constant density)
plt.plot(x_max_arr, x_max_result, linewidth = 1)
plt.scatter(x_max_arr, x_max_result, s = 12)
plt.xlabel("max wavenumber [$μm^{−1}$]", fontsize = 12)
plt.ylabel("log precision", fontsize = 12)
plt.xticks(fontsize = 12) 
plt.yticks(fontsize = 12) 
plt.yscale("log")
plt.savefig("Figures/hw_two_figures/sample_bii.png", bbox_inches='tight')
plt.show()

#precision improves as density increases 
plt.plot(density_arr, density_result, linewidth = 1)
plt.scatter(density_arr, density_result, s = 12)
plt.xlabel("density of wavenumber", fontsize = 12)
plt.ylabel("log precision", fontsize = 12)
plt.xticks(fontsize = 12) 
plt.yticks(fontsize = 12) 
plt.yscale("log")
plt.savefig("Figures/hw_two_figures/sample_biii.png", bbox_inches='tight')
plt.show()

"""
(c) Deliberately add a very small error to your calculation by making it overshoot the
right answer (just artificially add a constant to the output). What is the signature of such
an error, as you plot your precision vs. the three parameters, given the known answer?
"""
error = 1e-6

x_min_result = []
x_max_result = []
density_result = []
for i in x_min_arr:
    calculation = integrate_Planck(temp, i, 1e5, int(1e7)) + error
    precision = np.abs(truth - calculation)/truth
    x_min_result.append(precision)
for i in x_max_arr:
    calculation = integrate_Planck(temp, 0, i, int(1e7)) + error
    precision = np.abs(truth - calculation)/truth
    x_max_result.append(precision)
for i in density_arr:
    calculation = integrate_Planck(temp, 0, 1e5, int(i)) + error
    precision = np.abs(truth - calculation)/truth
    density_result.append(precision)

plt.plot(x_min_arr, x_min_result, linewidth = 1)
plt.scatter(x_min_arr, x_min_result, s = 12)
plt.xlabel("min wavenumber [$μm^{−1}$]", fontsize = 12)
plt.ylabel("log precision", fontsize = 12)
plt.xticks(fontsize = 12) 
plt.yticks(fontsize = 12) 
plt.yscale("log")
plt.savefig("Figures/hw_two_figures/sample_ci.png", bbox_inches='tight')
plt.show()

plt.plot(x_max_arr, x_max_result, linewidth = 1)
plt.scatter(x_max_arr, x_max_result, s = 12)
plt.xlabel("max wavenumber [$μm^{−1}$]", fontsize = 12)
plt.ylabel("log precision", fontsize = 12)
plt.xticks(fontsize = 12) 
plt.yticks(fontsize = 12) 
plt.yscale("log")
plt.savefig("Figures/hw_two_figures/sample_cii.png", bbox_inches='tight')
plt.show()

plt.plot(density_arr, density_result, linewidth = 1)
plt.scatter(density_arr, density_result, s = 12)
plt.xlabel("density of wavenumber", fontsize = 12)
plt.ylabel("log precision", fontsize = 12)
plt.xticks(fontsize = 12) 
plt.yticks(fontsize = 12) 
plt.yscale("log")
plt.savefig("Figures/hw_two_figures/sample_ciii.png", bbox_inches='tight')
plt.show()

"""
(d) Imagine you did not have the analytic answer to guide you in part (b). What tests can
you do to determine the precision of your integration, and to demonstrate that the values
for your three parameters are sufficient to achieve that precision? Show the results of
such a test.
"""
x_min_result = []
x_max_result = []
density_result = []
for i in x_min_arr:
    calculation = integrate_Planck(temp, i, 1e5, int(1e7))
    x_min_result.append(calculation)
for i in x_max_arr:
    calculation = integrate_Planck(temp, 0, i, int(1e7))
    x_max_result.append(calculation)
for i in density_arr:
    calculation = integrate_Planck(temp, 0, 1e5, int(i))
    density_result.append(calculation)

x_min_con = np.isclose(x_min_result, np.roll(x_min_result, shift=1), atol=1e-6)
x_max_con = np.isclose(x_max_result, np.roll(x_max_result, shift=1), atol=1e-6)
density_con = np.isclose(density_result, np.roll(density_result, shift=1), atol=1e-6)

precision_arr = []
x_min_index =  []
x_max_index =  []
density_index =  []
for i in range(0,len(x_min_con)-1):
    if (x_min_con[i] == True and x_min_con[i+1] == True):
        precision_arr.append(x_min_result[i])
        x_min_index.append(i)
        x_max_index.append(i)
        density_index.append(i)
    if (x_max_con[i] == True and x_max_con[i+1] == True):
        precision_arr.append(x_max_result[i])
        x_min_index.append(i)
        x_max_index.append(i)
        density_index.append(i)
    if (density_con[i] == True and density_con[i+1] == True):
        precision_arr.append(density_result[i])
        x_min_index.append(i)
        x_max_index.append(i)
        density_index.append(i)

precision_con = np.isclose(precision_arr, np.roll(precision_arr, shift=1), atol=1e-6)

final_precision_arr = []
final_x_min =  []
final_x_max =  []
final_density =  []
for i in range(1,len(precision_con)-1):
    if (precision_con[i] == True and precision_con[i+1] == True):
        final_precision_arr.append(precision_arr[i])
        final_x_min.append(x_min_arr[int(x_min_index[i])])
        final_x_max.append(x_max_arr[int(x_max_index[i])])
        final_density.append(density_arr[int(density_index[i])])

plt.plot(range(0,len(precision_arr)), precision_arr, linewidth = 1)
plt.scatter(range(0,len(precision_arr)), precision_arr, s = 12)
plt.scatter(0,integrate_Planck(temp, int(final_x_min[0]), int(final_x_max[0]), int(final_density[0])))
plt.xlabel("count", fontsize = 12)
plt.ylabel("intensity [erg $cm^{-2}$ $s^{-1}$ $st^{-1}$ $Hz^{-1}$]", fontsize = 12)
plt.xticks(fontsize = 12) 
plt.yticks(fontsize = 12) 
plt.yscale("log")
plt.savefig("Figures/hw_two_figures/sample_d.png", bbox_inches='tight')
plt.show()