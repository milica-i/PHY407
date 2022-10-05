#This code ___
# Authors: 

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
#2. ???

#the function we are considering for this question
def f(x, c):
    return 1 - np.e**(-c*x)

#Constants, initial guess, convergence threshold

c = 2 #for part a
x = 0.5 #initial guess
dx = 1 #initial distance
threshold = 1e-6
x_list = [x]


while dx > threshold:
    x_list.append(f(x_list[-1], c))
    dx = np.abs(x_list[-1] - x_list[-2])
print('The solution converges to', x_list[-1])




#Pseudocode for 3a - 6.10b
#1. Initialize c_array,
#2. ???

c_array = np.arange(0, 3, 0.01)
x_sols = np.zeros(len(c_array))

for i in range(len(c_array)):
    x = 0.5 #initial guess
    dx = 1 #initial distance
    x_list = [x]    
    while dx > threshold:
        x_list.append(f(x_list[-1], c_array[i]))
        dx = np.abs(x_list[-1] - x_list[-2]) 
    x_sols[i] = x_list[-1]
    
#print(x_sols)

plt.figure()
plt.plot(c_array, x_sols)
plt.xlabel('c')
plt.ylabel('x')
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
#Modify previous code to use overrelaxation 

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
    
print('The solution converges to', x_list[-1], 'after', counter, 'iterations.')



#Pseudocode for 3b - 6.11d
#Modify previous code to use overrelaxation 
def g(x):
    return 1 - np.e**(1-x**2)

x = 1 #initial guess
dx = 1 #initial distance
threshold = 1e-6
x_list = [x]
counter = 0


while dx > threshold:
    counter+= 1
    x_list.append(g(x_list[-1]))
    dx = np.abs(x_list[-1] - x_list[-2])
print('The solution converges to', x_list[-1], counter, 'iterations')

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
    print(counter, x_list[-1])
    
print('The solution converges to', x_list[-1], 'after', counter, 'iterations.')


#Question 3c
print('Question 3c')