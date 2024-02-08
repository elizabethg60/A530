from A530_package import box_integrator, intensity_integral
from A530_package import special_exp_integral, integrate_special_exp
import matplotlib.pyplot as plt
from scipy import integrate
import scipy.special as sc
import numpy as np

# 4. Emergent Intensity Integrator

"""
(a) Illustrate your answer to Problem 2 (from last week) by numerically integrating the
radial emergent intensity from a quadratic source function. That is, use
plots to show that the Eddington-Barbier approximation holds (i.e. its prediction for the value of Sν (τν = 1) is valid) 
for a linear source function, and that quadratic terms have the effects you derived in Problem 1.
"""
#hardwire values for a0 and a1
a0 = a1 = 1
#range of a2 values to test EB approximation
a2 = np.linspace(-1,1, 20)

#integral parameters (confirmed in 4b)
x_min = 10**(-7)
x_max = 50
step = 10**(-7)
x_array = np.arange(x_min, x_max, step)

#linear form 
a0_coeff = integrate.simpson(intensity_integral(0, x_array), x_array)
a1_coeff = integrate.simpson(intensity_integral(1, x_array), x_array)
output_linear = a0*a0_coeff+a1*a1_coeff
#quadratic form 
a2_coeff = integrate.simpson(intensity_integral(2, x_array), x_array)
output_quadratic = a0*a0_coeff+a1*a1_coeff+a2*a2_coeff
output_quadratic_zero = a0*a0_coeff+a1*a1_coeff+0*a2_coeff
source_quadratic = a0*a0_coeff+a1*a1_coeff+a2
#plot
plt.plot(a2, output_quadratic, linewidth = 2, color = 'k')
plt.plot(a2, source_quadratic, linewidth = 2, color = 'b')
plt.scatter(a2, output_quadratic, s = 15, color = 'k', label = "EB, Quadratic")
plt.scatter(a2, source_quadratic, s = 15, color = 'b', label = "Source, Quadratic")
plt.scatter(0, output_linear, s = 80, marker = "*", color = 'r', label = "Linear", zorder=5)
plt.xlabel("$a_{2}$", fontsize = 12)
plt.ylabel("[erg $cm^{-2}$ $s^{-1}$ $st^{-1}$ $Hz^{-1}$]", fontsize = 12)
plt.xticks(fontsize = 12) 
plt.yticks(fontsize = 12)
plt.legend()
plt.savefig("Figures/hw_three_figures/sample_4a.pdf", bbox_inches='tight')
plt.show()

#confirm: negative values of a2 results in EB underestimate whereas positive values of a2 results in EB overestimate and at a2 = 0 perfect match

"""
(b) Explain (and perhaps show) how you determined that your sampling and limits of
integration in part (a) were appropriate and how you determined your precision. (When
quantifying precision, what is the appropriate metric?)
"""
#hardwire value for a2
a2 = 2 

#precision for linear and quadratic functions with optimized integral parameters (above)
truth_linear = a0+a1
truth_quadratic = a0+a1+2*a2
output_quadratic = a0*a0_coeff+a1*a1_coeff+a2*a2_coeff
precision_linear = np.abs(truth_linear - output_linear)/truth_linear
precision_quadratic = np.abs(truth_quadratic - output_quadratic)/truth_quadratic

print("truth: {} vs calculation: {} gives precision of {}".format(truth_linear, output_linear, precision_linear))
#truth: 2 vs calculation: 2.0 gives precision of 0.0
#truth: 2 vs calculation: 1.9999998999999513 gives precision of 5.00000243430776e-08
print("truth: {} vs calculation: {} gives precision of {}".format(truth_quadratic, output_quadratic, precision_quadratic))
#truth: 6 vs calculation: 5.999999999999998 gives precision of 2.9605947323337506e-16
#truth: 6 vs calculation: 5.9999998999998425 gives precision of 1.66666929146686e-08

#relative prescision convergence
# x_min_arr = np.logspace(-20, -10, 20)
# x_max_arr = np.logspace(1, 2, 20)
# step_arr = np.logspace(-6, -5, 20)

x_min_arr = np.logspace(-7, -6, 20)
x_max_arr = np.logspace(1, 2, 20)
step_arr = np.logspace(-7, -6, 20)
 
x_min_result_quadratic = []
x_max_result_quadratic = []
step_result_quadratic = []
for i in x_min_arr:
    x_array = np.power(10, np.arange(np.log10(i), np.log10(x_max), step))
    a0_coeff = integrate.simpson(intensity_integral(0, x_array), x_array)
    a1_coeff = integrate.simpson(intensity_integral(1, x_array), x_array)
    a2_coeff = integrate.simpson(intensity_integral(2, x_array), x_array)
    output_quadratic = a0*a0_coeff+a1*a1_coeff+a2*a2_coeff
    precision_quadratic = np.abs(truth_quadratic - output_quadratic)/truth_quadratic
    x_min_result_quadratic.append(precision_quadratic)
for i in x_max_arr:
    x_array = np.power(10, np.arange(np.log10(x_min), np.log10(i), step))
    a0_coeff = integrate.simpson(intensity_integral(0, x_array), x_array)
    a1_coeff = integrate.simpson(intensity_integral(1, x_array), x_array)
    a2_coeff = integrate.simpson(intensity_integral(2, x_array), x_array)
    output_quadratic = a0*a0_coeff+a1*a1_coeff+a2*a2_coeff
    precision_quadratic = np.abs(truth_quadratic - output_quadratic)/truth_quadratic
    x_max_result_quadratic.append(precision_quadratic)
for i in step_arr:
    x_array = np.power(10, np.arange(np.log10(x_min), np.log10(x_max), i))
    a0_coeff = integrate.simpson(intensity_integral(0, x_array), x_array)
    a1_coeff = integrate.simpson(intensity_integral(1, x_array), x_array)
    a2_coeff = integrate.simpson(intensity_integral(2, x_array), x_array)
    output_quadratic = a0*a0_coeff+a1*a1_coeff+a2*a2_coeff
    precision_quadratic = np.abs(truth_quadratic - output_quadratic)/truth_quadratic
    step_result_quadratic.append(precision_quadratic)

plt.plot(x_min_arr, x_min_result_quadratic, linewidth = 2, color = 'k')
plt.scatter(x_min_arr, x_min_result_quadratic, s = 15, color = 'k')
plt.xlabel("log min bound", fontsize = 12)
plt.ylabel("log relative precision", fontsize = 12)
plt.xticks(fontsize = 12) 
plt.yticks(fontsize = 12) 
plt.yscale("log")
plt.xscale("log")
plt.legend()
plt.savefig("Figures/hw_three_figures/sample_4i.pdf", bbox_inches='tight')
plt.show()

plt.plot(x_max_arr, x_max_result_quadratic, linewidth = 2, color = 'k')
plt.scatter(x_max_arr, x_max_result_quadratic, s = 15, color = 'k')
plt.xlabel("log max bound", fontsize = 12)
plt.ylabel("log relative precision", fontsize = 12)
plt.xticks(fontsize = 12) 
plt.yticks(fontsize = 12) 
plt.yscale("log")
plt.xscale("log")
plt.legend()
plt.savefig("Figures/hw_three_figures/sample_4bii.pdf", bbox_inches='tight')
plt.show()

plt.plot(step_arr, step_result_quadratic, linewidth = 2, color = 'k')
plt.scatter(step_arr, step_result_quadratic, s = 15, color = 'k')
plt.xlabel("log step width", fontsize = 12)
plt.ylabel("log precision", fontsize = 12)
plt.xticks(fontsize = 12) 
plt.yticks(fontsize = 12) 
plt.yscale("log")
plt.xscale("log")
plt.legend()
plt.savefig("Figures/hw_three_figures/sample_4biii.pdf", bbox_inches='tight')
plt.show()


## 5. Exponential Integrals

"""
a) Write or familiarize yourself with Python's scipy.special.expn or the equivalent
exponential integral functions in whatever language you are using
"""
#using scipy.special.expn 
#inputs: n - arrray of non negative integers and x array of real non negative integers
#outputs: values of generalized exponential integral 

x = np.array([1, 2, 3, 4])
print("E(n = 0) = {}".format(sc.expn(0, x)))
print("E(n = 1) = {}".format(sc.expn(1, x)))

"""
b) Numerically integrate the area under the first three exponential functions (that is,
numerically integrate En(x)dx for n = 1, 2, 3). How close can you get to the analytic result? What limits your precision?
"""

#analytic results (derived via mathematica - shared in report)
n = [1,2,3]
truth1 = 1/n[0]
truth2 = 1/n[1]
truth3 = 1/n[2]
#numerically integrate 
x_min = 10**(-7) 
x_max = 50
step = 10**(-7) 
x_array = np.power(10, np.arange(np.log10(x_min), np.log10(x_max), step))
expn1 = sc.expn(1, x_array)
expn2 = sc.expn(2, x_array)
expn3 = sc.expn(3, x_array)
output1 = integrate.simpson(expn1, x_array)
output2 = integrate.simpson(expn2, x_array)
output3 = integrate.simpson(expn3, x_array)
#precision for n = 1,2,3 with optimized integral parameters
precision1 = np.abs(truth1 - output1)/truth1
precision2 = np.abs(truth2 - output2)/truth2
precision3 = np.abs(truth3 - output3)/truth3

print("truth: {} vs calculation: {} gives precision of {}".format(truth1, output1, precision1))
#truth: 1.0 vs calculation: 0.9999999999999962 gives precision of 3.774758283725532e-15
#truth: 1.0 vs calculation: 0.9999983459119953 gives precision of 1.6540880046767015e-06
print("truth: {} vs calculation: {} gives precision of {}".format(truth2, output2, precision2))
#truth: 0.5 vs calculation: 0.4999999999999999 gives precision of 2.220446049250313e-16
#truth: 0.5 vs calculation: 0.4999999000000842 gives precision of 1.999998315849183e-07
print("truth: {} vs calculation: {} gives precision of {}".format(truth3, output3, precision3))
#truth: 0.3333333333333333 vs calculation: 0.33333333333333337 gives precision of 1.6653345369377348e-16
#truth: 0.3333333333333333 vs calculation: 0.33333328333333756 gives precision of 1.4999998726450414e-07

"""
c) Explain how you determined that your sampling and upper limit in part (b) are appro-
priate and how you determined the precision of your integration. For a box or trapezoid
integrator, you may need to experiment with appropriate samplings and limits to achieve
an precise result. What’s the best precision you can get in a reasonable computation time
(less than a second)? What upper limit and spacing of points gives you that precision?
How few points can you get away with, practically speaking? (If you use a black-box
integrator, describe how it does these things and how you know it’s working.)
"""

#given n = 1 had the worse precision with the parameters above, i will optimize parameters using n =1 

#relative prescision convergence
# x_min_arr = np.logspace(-20, -10, 20)
# x_max_arr = np.logspace(1, 2, 20)
# step_arr = np.logspace(-6, -5, 20)
 
x_min_arr = np.logspace(-7, -6, 20)
x_max_arr = np.logspace(1, 2, 20)
step_arr = np.logspace(-7, -6, 20)

x_min_result = []
x_max_result = []
step_result = []
for i in x_min_arr:
    x_array = np.power(10, np.arange(np.log10(i), np.log10(x_max), step))
    expn1 = sc.expn(1, x_array)
    output1 = integrate.simpson(expn1, x_array)
    precision1 = np.abs(truth1 - output1)/truth1
    x_min_result.append(precision1)
for i in x_max_arr:
    x_array = np.power(10, np.arange(np.log10(x_min), np.log10(i), step))
    expn1 = sc.expn(1, x_array)
    output1 = integrate.simpson(expn1, x_array)
    precision1 = np.abs(truth1 - output1)/truth1
    x_max_result.append(precision1)
for i in step_arr:
    x_array = np.power(10, np.arange(np.log10(x_min), np.log10(x_max), i))
    expn1 = sc.expn(1, x_array)
    output1 = integrate.simpson(expn1, x_array)
    precision1 = np.abs(truth1 - output1)/truth1
    step_result.append(precision1)

plt.plot(x_min_arr, x_min_result, linewidth = 2, color = 'k')
plt.scatter(x_min_arr, x_min_result, s = 15, color = 'k')
plt.xlabel("log min bound", fontsize = 12)
plt.ylabel("log relative precision", fontsize = 12)
plt.xticks(fontsize = 12) 
plt.yticks(fontsize = 12) 
plt.yscale("log")
plt.xscale("log")
plt.legend()
plt.savefig("Figures/hw_three_figures/sample_5bi.pdf", bbox_inches='tight')
plt.show()

plt.plot(x_max_arr, x_max_result, linewidth = 2, color = 'k')
plt.scatter(x_max_arr, x_max_result, s = 15, color = 'k')
plt.xlabel("log max bound", fontsize = 12)
plt.ylabel("log relative precision", fontsize = 12)
plt.xticks(fontsize = 12) 
plt.yticks(fontsize = 12) 
plt.yscale("log")
plt.xscale("log")
plt.legend()
plt.savefig("Figures/hw_three_figures/sample_5bii.pdf", bbox_inches='tight')
plt.show()

plt.plot(step_arr, step_result, linewidth = 2, color = 'k')
plt.scatter(step_arr, step_result, s = 15, color = 'k')
plt.xlabel("log step width", fontsize = 12)
plt.ylabel("log precision", fontsize = 12)
plt.xticks(fontsize = 12) 
plt.yticks(fontsize = 12) 
plt.yscale("log")
plt.xscale("log")
plt.legend()
plt.savefig("Figures/hw_three_figures/sample_5biii.pdf", bbox_inches='tight')
plt.show()

#numerically integrate 
x_min = 10**(-20) 
x_max = 50
step = 10**(-5) 
x_array = np.power(10, np.arange(np.log10(x_min), np.log10(x_max), step))
expn1 = sc.expn(1, x_array)
expn2 = sc.expn(2, x_array)
expn3 = sc.expn(3, x_array)
output1 = integrate.simpson(expn1, x_array)
output2 = integrate.simpson(expn2, x_array)
output3 = integrate.simpson(expn3, x_array)
#precision for n = 1,2,3 with optimized integral parameters
precision1 = np.abs(truth1 - output1)/truth1
precision2 = np.abs(truth2 - output2)/truth2
precision3 = np.abs(truth3 - output3)/truth3

print("truth: {} vs calculation: {} gives precision of {}".format(truth1, output1, precision1))
#truth: 1.0 vs calculation: 1.0 gives precision of 0.0
print("truth: {} vs calculation: {} gives precision of {}".format(truth2, output2, precision2))
#truth: 0.5 vs calculation: 0.5 gives precision of 0.0
print("truth: {} vs calculation: {} gives precision of {}".format(truth3, output3, precision3))
#truth: 0.3333333333333333 vs calculation: 0.33333333333333326 gives precision of 1.6653345369377348e-16