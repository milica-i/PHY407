# The following code evaluates a given integral using the Trapezoidal rule 
# and Simpson's method for integration and compares and contrasts the results.
# Author: Milica Ivetic

import numpy as np
from time import time
import trapezoidal_simpson as trap_sim

# Question 2b

# This is the function being considered in this question

def f(x):
    return 4 / (1+x**2)

# values from using trapezoidal rule and simpson's rule for N=4

trap_b = trap_sim.trapezoidal(f, 4, 0, 1)
simpson_b = trap_sim.simpson(f, 4, 0, 1)

print('Question 2b')
print('The value obtained for N=4 using Trapezoidal rule is', trap_b)
print('The value obtained for N=4 using Simpson\'s rule is', simpson_b)
print('The exact value is', np.pi)


# Question 2c

#Considering the Trapezoidal method:

#The calculation for the error (difference between exact value and calculated value)

print('Question 2c')

print('The error for trapezoidal rule for N=10000 is:',np.pi-trap_sim.trapezoidal(f, 10000, 0, 1))

#Setting up timer to time how long this computation takes
start_trap = time()

trap_sim.trapezoidal(f, 10000, 0, 1)

end_trap = time()

diff_trap = end_trap - start_trap

print('The calculation for the trapezoidal rule for N=10000 takes', diff_trap, 'seconds to compute.')



#Considering Simpson's method:

print('The error for Simpson\'s rule for N=14 is:',np.pi-trap_sim.simpson(f, 14, 0, 1))


#Setting up timer to time how long this computation takes
start_simpson = time()

#integral calculation here
trap_sim.simpson(f, 14, 0, 1)

end_simpson = time()

diff_simpson = end_simpson - start_simpson

print('The calculation for Simpson\'s rule for N=14 takes', diff_simpson, 'seconds to compute.')



# Question 2d

print('Question 2d')

#Using "practical estimation of errors" method from section 5.2.1, p. 153 of the textbook

#the following is given in the question 
N1 = 16
N2 = 32

# Value of integral using N1
I1 = trap_sim.trapezoidal(f, N1, 0, 1)

#Value of integral using N2
I2 = trap_sim.trapezoidal(f, N2, 0, 1)


#Using equation 5.28 from the textbook

error_2 = (1/3)*(I2 - I1)

print('The practical error estimation for Trapezoidal rule with N2 = 32 is', error_2)


