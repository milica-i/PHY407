# This code considers different methods for solving linear systems:
# Gaussian Elimination, Partial Pivoting, and LU Decomposition. 
# We compare how long each method takes to run for various array sizes, and 
# we analyze their errors.

#Author: Milica Ivetic


import numpy as np
import matplotlib.pyplot as plt
from SolveLinear import PartialPivot, GaussElim
from time import time


plt.rcParams['axes.titlesize'] = 13
plt.rcParams['axes.labelsize'] = 13
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12

#Question 1a
print('Question 1a')
print('Checking if PartialPivot code returns answer 6.16 to Equation 6.2 of Newman.')

Atest = np.array([[2, 1, 4, 1], [3,4,-1, -1],[1,-4,1, 5],[2,-2,1,3]], float)
Vtest = np.array([-4, 3, 9,7], float)
print('Partial Pivot returns', PartialPivot(Atest, Vtest),'. Which matches the answer in Newman.')

#Question 1b
print('Question 1b - see plots')

#Pseudocode for 1b
#1. Initialize N array. For each of a range of values of N, creates a random matrix A and a random array v. Also initialize timings and error arrays
#2. Find the solution x for the same A and v for the three different methods. 
    #Also measure the time it takes to solve for x using each method using time package
#3. For each case check the answer by comparing vsol = np.dot(A,x) vs. v by
    #calculating the mean of abs differences between the arrays
#4. plot N vs errors and N vs timings for each method


#1. Initialize N array. For each of a range of values of N, creates a random matrix A and a random array v. Also initialize timings and error arrays.

N_array = np.linspace(5, 500, 50)

gauss_timings = np.zeros(len(N_array))
pivot_timings = np.zeros(len(N_array))
LU_timings = np.zeros(len(N_array))

gauss_errors = np.zeros(len(N_array))
pivot_errors = np.zeros(len(N_array))
LU_errors = np.zeros(len(N_array))

#Steps 2 and 3
for i in range(len(N_array)):
    v = np.random.rand(int(N_array[i]))
    A = np.random.rand(int(N_array[i]), int(N_array[i]))
    
    #Using GaussElim
    
    gauss_start = time()
    
    gauss_sol = GaussElim(A, v) #perform calculation
    
    gauss_end = time()
    gauss_diff = gauss_end - gauss_start
    
    gauss_timings[i] = gauss_diff #add timing to time array
    
    #Calculate error and add it to error array
    
    gauss_vsol = np.dot(A, gauss_sol)
    gauss_err = np.mean(np.abs(v - gauss_vsol))
    
    gauss_errors[i] = gauss_err
    
    
    
    #Using PartialPivot
    
    pivot_start = time()
    
    pivot_sol = PartialPivot(A, v) #perform calculation
    
    pivot_end = time()
    pivot_diff = pivot_end - pivot_start
    
    pivot_timings[i] = pivot_diff  #add timing to time array
    
    #Calculate error and add it to error array
    
    pivot_vsol = np.dot(A, pivot_sol)
    pivot_err = np.mean(np.abs(v - pivot_vsol))
    
    pivot_errors[i] = pivot_err    
    
    
    #Using LU decomposition, which is done using np.linalg.solve
    
    LU_start = time()
    
    LU_sol = np.linalg.solve(A, v) # perform calculation 
                      
    LU_end = time()
    LU_diff = LU_end - LU_start
    
    LU_timings[i] = LU_diff   #add timing to time array
    
    #Calculate error and add it to error array
    
    LU_vsol = np.dot(A, LU_sol)
    LU_err = np.mean(np.abs(v - LU_vsol))
    
    LU_errors[i] = LU_err
    
    
#4. plot N vs errors and N vs timings for each method

plt.figure()
plt.plot(N_array, (gauss_timings), label = 'GaussElim')
plt.plot(N_array, (pivot_timings), label = 'PartialPivot')
plt.plot(N_array, (LU_timings), label = 'LU Decomposition')
plt.yscale('log')
plt.legend()
plt.xlabel('N')
plt.ylabel('Timing')
plt.title('Times to run each method as a function of N')
plt.show()

plt.figure()
plt.plot(N_array, (gauss_errors), label = 'GaussElim')
plt.plot(N_array, (pivot_errors), label = 'PartialPivot')
plt.plot(N_array, (LU_errors), label = 'LU Decomposition')
plt.yscale('log')
plt.legend()
plt.xlabel('N')
plt.ylabel('Error')
plt.title('Errors from each method as a function of N')
plt.show()
    