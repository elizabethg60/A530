import numpy as np
import pandas as pd
from A530_package import Pe_equation

"""
Write a program to implement Equation (9.8) from D. F. Gray (3rd edition) for parts
of the Sun. Consider the all elements for which you have partition functions, ionizations,
and abundances (but you may ignore second ionizations if you like).
You will need a table of abundances and atomic weights; use the solar abundances file in
the Solar Atmospheric Data folder on Canvas. This table comes from Table 16.3 of Gray
(3rd ed.).
Since you will want to use this program multiple times in later problems, you ultimately
should put the computational burden into a subroutine that takes the abundances, T ,
and Pg as input and returns Pe. Gray has a good discussion on initial estimates for Pe.
Use Table 9.2 in Gray to check your code.
Please be verbose. For instance, you should describe how your program chooses an initial
estimate of Pe, based on the input values of (T , log Pg), how your program chooses
subsequent estimates of Pe from prior iterations, and the criterion that determines when
to stop iterating. If you use a black box, you still need to describe how it does these
things and how you checked that you understand it and that it’s working properly.
"""


temperature = 4000
Pe = 10
Pg = 5

print(Pe_equation(Pg, temperature, Pe))

