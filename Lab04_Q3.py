#This code explores different methods for solving non-linear systems:
# Relaxation, Overrelaxation, and Binary Search. Questions from Newman (2012)

# Author: Milica Ivetic

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['axes.titlesize'] = 13
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12


#Question 3a
print('Question 3a')

#Pseudocode for 3a - 6.10a
#1. Define function, constants, initial guess, and convergence threshold, x_list
#2. Perform relaxation calculation using method described in Newman

#1.
#the function we are considering for this question
def f(x, c):
    return 1 - np.e**(-c*x)

#Constants, initial guess, convergence threshold

c = 2 #for part a
x = 0.5 #initial guess
dx = 1 #initial distance
threshold = 1e-6
x_list = [x]

#2. Relaxation method as described in Newman
while dx > threshold:
    x_list.append(f(x_list[-1], c)) #update list of solutions
    dx = np.abs(x_list[-1] - x_list[-2]) #update distance
print('The solution converges to', x_list[-1])




#Pseudocode for 3a - 6.10b
#1. Initialize c_array, and xsols array
#2. For each c value, perform relaxation calculation, save results to array
#3. Plot solution as a function of c

#1. Initialize c_array, and xsols array
c_array = np.arange(0, 3, 0.01)
x_sols = np.zeros(len(c_array))


#2. For each c value, perform relaxation calculation, save results to array
for i in range(len(c_array)):
    x = 0.5 #initial guess
    dx = 1 #initial distance
    x_list = [x]    
    #same code as before, just varying c value
    while dx > threshold:
        x_list.append(f(x_list[-1], c_array[i]))
        dx = np.abs(x_list[-1] - x_list[-2]) 
    x_sols[i] = x_list[-1]
    

#3. Plot solution as a function of c
plt.figure()
plt.plot(c_array, x_sols)
plt.xlabel('c')
plt.ylabel('x')
plt.title('Solution as a function of c')
plt.show()



#Question 3b
print('Question 3b')

#Pseudocode for 3b - 6.11b
#Same code as in 6.10, just modified to print the number of iterations for c=2

c = 2 
x = 0.5 #initial guess
dx = 1 #initial distance
threshold = 1e-6
x_list = [x]
counter =0 #adding iteration counter

while dx > threshold:
    counter += 1
    print(counter)
    x_list.append(f(x_list[-1], c))
    dx = np.abs(x_list[-1] - x_list[-2])
print('The solution converges to', x_list[-1], 'after', counter, 'iterations.')


#Pseudocode for 3b - 6.11c
#Modify previous code to use overrelaxation method

c = 2 
x = 0.5 #initial guess
dx = 1 #initial distance
threshold = 1e-6
x_list = [x]
counter =0 #adding iteration counter
w = 0.5 #starting point for this parameter

while dx > threshold:
    counter += 1
    x_list.append((1+w)*f(x_list[-1], c) - w*x_list[-1])
    dx = np.abs(x_list[-1] - x_list[-2])
    print(counter, x_list[-1])
    
print('The solution converges to', x_list[-1], 'after', counter, 'iterations using overrelaxation method.')


#IGNORE 6.11d
#Pseudocode for 3b - 6.11d
#Modify previous code to use overrelaxation 
def g(x):
    return 1 - np.e**(1-x**2)

x = 0.5 #initial guess
dx = 1 #initial distance
threshold = 1e-6
x_list = [x]
counter = 0


while dx > threshold:
    counter+= 1
    x_list.append(g(x_list[-1]))
    dx = np.abs(x_list[-1] - x_list[-2])
#print('The solution converges to', x_list[-1], counter, 'iterations')
#print(counter)

x = 0.5 #initial guess
dx = 1 #initial distance
threshold = 1e-6
x_list = [x]
counter =0 #adding iteration counter
w = -0.5 #starting point for this parameter

while dx > threshold:
    counter += 1
    x_list.append((1+w)*g(x_list[-1]) - w*x_list[-1])
    dx = np.abs(x_list[-1] - x_list[-2])
    
#print('The solution converges to', x_list[-1], 'after', counter, 'iterations \
 #for x = e^(1-x^2) ')
#print(counter)


#Question 3c
print('Question 3c')

#Pseudocode for 3c - 6.13b
#1. Define function we are considering for this question and accuracy threshold
#2. Perform binary search method as outlined in Newman:

#given an initial pair of points, check that f(x) have opposite signs
# calculate the midpoint and eval f(mid)
#if f(mid) has the same sign as f(x1), then x1 = mid, otherwise, x2 = mid
#if abs(x1 - x2) > accuracy, repeat from step 2, otherwise, calculate new midpoint
# once more and this is the final estimate of the position of the root

#3. Calculate value for displacement constant


#1. Define function we are considering for this question, constants, and accuracy threshold

acc = 1e-6 #accuracy threshold

#constants:
c = 3e8 #m/s, speed of light
h_const = 6.62607015e-34 #m^2kg/s, Planck's constant
k = 1.380649e-23 #m^2 kg s^-2 K^-1, Boltzmann's constant

def h(x):
    return (-5)*(np.e**(-x)) + 5 -x #based on result from part 6.13a

#2. Perform binary search method as outlined in Newman:

#given an initial pair of points, check that h(x1) and h(x2) have opposite signs
#note that we do NOT want the obvious root at x = 0

x1 = 1
x2 = 6

print('x1 =', x1, 'x2 =', x2)
print('h(x1) = ', h(x1), 'h(x2) = ', h(x2))
print('The function values at the two initial guess points have opposite signs.')

#in a while loop:
# calculate the midpoint and eval h(mid)
#if h(mid) has the same sign as h(x1), then x1 = mid, otherwise, x2 = mid

#if abs(x1 - x2) > accuracy, repeat from step 2, otherwise, calculate new midpoint once more
while np.abs(x1 - x2) > acc:
    mid = 0.5*(x1 + x2) #midpoint
    func_mid = h(mid) #function at midpoint
   #check signs of function at midpoint and at x1
    if np.sign(func_mid) == np.sign(h(x1)):
        x1 = mid
    else:
        x2 = mid
        
final_midpoint = 0.5*(x1 + x2)
print('The final estimation is', final_midpoint)

#3. Calculate value for displacement constant
b = h_const*c/(k*final_midpoint)
    
    
print('The displacement constant is', b)


#Pseudocode for 3c - 6.13c

#1. Define lambda 
#2. Using lambda = b/T, solve for T


#1. Define lambda

l = 5.02e-7 #m , wavelength


#2. #2. Using lambda = b/T, solve for T

T = b/l

print('The estimated surface temperature of the Sun is', T, 'K')