from A530_package import box_integrator, intensity_integral
from A530_package import special_exp_integral, integrate_special_exp
import matplotlib.pyplot as plt
from scipy import integrate
import scipy.special as sc
import numpy as np

## 4. Emergent Intensity Integrator

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
x_min = 10**(-20)
x_max = 50
step = 10**(-5)
x_array = np.arange(x_min, x_max, step)

#linear form 
a0_coeff = integrate.simpson(intensity_integral(0, x_array), x_array)
a1_coeff = integrate.simpson(intensity_integral(1, x_array), x_array)
output_linear = a0*a0_coeff+a1*a1_coeff
#quadratic form 
a2_coeff = integrate.simpson(intensity_integral(2, x_array), x_array)
output_quadratic = a0*a0_coeff+a1*a1_coeff+a2*a2_coeff
output_quadratic_zero = a0*a0_coeff+a1*a1_coeff+0*a2_coeff
#plot
plt.plot(a2, [output_linear]*len(a2), linewidth = 2, color = 'r')
plt.plot(a2, output_quadratic, linewidth = 2, color = 'k')
plt.scatter(a2, [output_linear]*len(a2), s = 15, color = 'r', label = "Linear")
plt.scatter(a2, output_quadratic, s = 15, color = 'k', label = "Quadratic")
plt.scatter(0, output_quadratic_zero, s = 15, color = 'k')
plt.xlabel("$a_{2}$", fontsize = 12)
plt.ylabel("emergent intensity [erg $cm^{-2}$ $s^{-1}$ $st^{-1}$ $Hz^{-1}$]", fontsize = 12)
plt.xticks(fontsize = 12) 
plt.yticks(fontsize = 12)
plt.legend()
plt.savefig("Figures/hw_three_figures/sample_4a.png", bbox_inches='tight')
plt.show()

#confirm: negative values of a2 results in EB overestimating whereas positive values of a2 results in EB underestimating and at a2 = 0 perfect match

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
print("truth: {} vs calculation: {} gives precision of {}".format(truth_quadratic, output_quadratic, precision_quadratic))
#truth: 6 vs calculation: 5.999999999999998 gives precision of 2.9605947323337506e-16

#relative prescision convergence
x_min_arr = np.logspace(-20, -10, 20)
x_max_arr = np.logspace(1, 2, 20)
step_arr = np.logspace(-6, -5, 20)
 
x_min_result_quadratic = []
x_max_result_quadratic = []
step_result_quadratic = []
for i in x_min_arr:
    x_array = np.arange(i, x_max, step)
    a0_coeff = integrate.simpson(intensity_integral(0, x_array), x_array)
    a1_coeff = integrate.simpson(intensity_integral(1, x_array), x_array)
    a2_coeff = integrate.simpson(intensity_integral(2, x_array), x_array)
    output_quadratic = a0*a0_coeff+a1*a1_coeff+a2*a2_coeff
    precision_quadratic = np.abs(truth_quadratic - output_quadratic)/truth_quadratic
    x_min_result_quadratic.append(precision_quadratic)
for i in x_max_arr:
    x_array = np.arange(x_min, i, step)
    a0_coeff = integrate.simpson(intensity_integral(0, x_array), x_array)
    a1_coeff = integrate.simpson(intensity_integral(1, x_array), x_array)
    a2_coeff = integrate.simpson(intensity_integral(2, x_array), x_array)
    output_quadratic = a0*a0_coeff+a1*a1_coeff+a2*a2_coeff
    precision_quadratic = np.abs(truth_quadratic - output_quadratic)/truth_quadratic
    x_max_result_quadratic.append(precision_quadratic)
for i in step_arr:
    x_array = np.arange(x_min, x_max, i)
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
plt.savefig("Figures/hw_three_figures/sample_4i.png", bbox_inches='tight')
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
plt.savefig("Figures/hw_three_figures/sample_4bii.png", bbox_inches='tight')
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
plt.savefig("Figures/hw_three_figures/sample_4biii.png", bbox_inches='tight')
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
#precision for linear and quadratic functions with optimized integral parameters (above)
truth_linear = a0+a1
truth_quadratic = a0+a1+2*a2
output_quadratic = a0*a0_coeff+a1*a1_coeff+a2*a2_coeff
precision_linear = np.abs(truth_linear - output_linear)/truth_linear
precision_quadratic = np.abs(truth_quadratic - output_quadratic)/truth_quadratic

print("truth: {} vs calculation: {} gives precision of {}".format(truth_linear, output_linear, precision_linear))
#truth: 2 vs calculation: 2.0 gives precision of 0.0
print("truth: {} vs calculation: {} gives precision of {}".format(truth_quadratic, output_quadratic, precision_quadratic))
#truth: 6 vs calculation: 5.999999999999998 gives precision of 2.9605947323337506e-16

#relative prescision convergence
x_min_arr = np.logspace(-20, -10, 20)
x_max_arr = np.logspace(1, 2, 20)
step_arr = np.logspace(-6, -5, 20)
 
x_min_result_quadratic = []
x_max_result_quadratic = []
step_result_quadratic = []
for i in x_min_arr:
    x_array = np.arange(i, x_max, step)
    a0_coeff = integrate.simpson(intensity_integral(0, x_array), x_array)
    a1_coeff = integrate.simpson(intensity_integral(1, x_array), x_array)
    a2_coeff = integrate.simpson(intensity_integral(2, x_array), x_array)
    output_quadratic = a0*a0_coeff+a1*a1_coeff+a2*a2_coeff
    precision_quadratic = np.abs(truth_quadratic - output_quadratic)/truth_quadratic
    x_min_result_quadratic.append(precision_quadratic)
for i in x_max_arr:
    x_array = np.arange(x_min, i, step)
    a0_coeff = integrate.simpson(intensity_integral(0, x_array), x_array)
    a1_coeff = integrate.simpson(intensity_integral(1, x_array), x_array)
    a2_coeff = integrate.simpson(intensity_integral(2, x_array), x_array)
    output_quadratic = a0*a0_coeff+a1*a1_coeff+a2*a2_coeff
    precision_quadratic = np.abs(truth_quadratic - output_quadratic)/truth_quadratic
    x_max_result_quadratic.append(precision_quadratic)
for i in step_arr:
    x_array = np.arange(x_min, x_max, i)
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
plt.savefig("Figures/hw_three_figures/sample_4i.png", bbox_inches='tight')
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
plt.savefig("Figures/hw_three_figures/sample_4bii.png", bbox_inches='tight')
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
plt.savefig("Figures/hw_three_figures/sample_4biii.png", bbox_inches='tight')
plt.show()















# # #using scipy.special.expn to get analytic result for n = 1,2,3 for x = 10**(-10) (which is lower bound of integral)
# # # x_min = 10**(-10)
# # # truth_E0 = sc.expn(0, x_min)
# # # truth_E1 = sc.expn(1, x_min)
# # # truth_E2 = sc.expn(2, x_min)








# # # x_min_arr = np.logspace(-12, -5, 50)
# # # x_max_arr = np.logspace(2, 3, 20)
# # # step_arr = np.logspace(-6, -5, 20)

# # # x_min_result = []
# # # x_max_result = []
# # # step_result = []

# # # x_min = 10**(-4)
# # # x_max = 100
# # # step = 0.00001

# # # # for i in x_min_arr:
# # # #     output = integrate_special_exp(0, i, x_max, step)
# # # #     precision = np.abs(truth_E0 - output)/truth_E0
# # # #     x_min_result.append(precision)
# # # for i in x_max_arr:
# # #     output = integrate_special_exp(0, x_min, i, step)
# # #     precision = np.abs(truth_E0 - output)/truth_E0
# # #     x_max_result.append(precision)
# # # # for i in step_arr:
# # # #     output = integrate_special_exp(0, x_min, x_max, i)
# # # #     precision = np.abs(truth_E0 - output)/truth_E0
# # # #     step_result.append(precision)

# # # # plt.plot(x_min_arr, x_min_result, linewidth = 1)
# # # # plt.xticks(fontsize = 12) 
# # # # plt.yticks(fontsize = 12) 
# # # # plt.yscale("log")
# # # # plt.xscale("log")
# # # # plt.show()

# # # plt.plot(x_max_arr, x_max_result, linewidth = 1)
# # # plt.xticks(fontsize = 12) 
# # # plt.yticks(fontsize = 12) 
# # # plt.yscale("log")
# # # plt.xscale("log")
# # # plt.show()

# # # # plt.plot(step_arr, step_result, linewidth = 1)
# # # # plt.xticks(fontsize = 12) 
# # # # plt.yticks(fontsize = 12) 
# # # # plt.yscale("log")
# # # # plt.xscale("log")
# # # # plt.show()




# # # # #print(integrate_special_exp(0, x_min, 100000, 0.1))









