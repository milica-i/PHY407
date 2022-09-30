# This code uses Gaussian quadrature to numerically calculate the period of the spring with the period
#given by eqn 7 in the lab manual, and see how it transitions from the classical to the relativistic case.
# Author: Milica Ivetic


#import packages and define constants

import numpy as np
import matplotlib.pyplot as plt
from gaussxw import gaussxwab

plt.rcParams['axes.titlesize'] = 13
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12

m = 1 #kg, mass
k = 12 #N/m, spring constant
c = 3e8 #m/s, speed of light
x0 = 0.01 #m, initial position we are considering for this problem, in metres 
          #to stay consistent with units

# Question 2a
print('Question 2a')

#Pseudocode for Q2a:
#1. import packages and define constants
#2. Calculate the classical value 2pi*sqrt(m/k)
#3. Define function g(x), given by Eqn. 7 of the lab
#4.Use Gaussian quadrature method to find the period for N=8
#4.1 Need to calculate the sample points and weights, then map them, will use gaussxwab
#This function is more convenient because it does the mapping for us 
#4.2 Perform the integration, output calculated period value
#5. Do the same for N=16
#6. Estimate error using Eqn 5.66 of the textbook AND by comparing with the classical value


#2. Calculate the classical value 2pi*sqrt(m/k)

classical_T = 2*np.pi*np.sqrt(m/k)
print('The classical value of the period is:', classical_T, 'seconds')

#3. Define function g(x), given by Eqn. 7 of the lab. Define f(x) = 1/g(x)

def g(x, x0):
    diff = x0**2 - x**2
    return c*((k*(diff)*(2*m*c**2+k*(diff)/2))/(2*(m*c**2+k*(diff)/2)**2))**(0.5)

def f(x, x0):
    return 1/g(x,x0)

#4.Use Gaussian quadrature method to find the period for N=8

#4.1 Need to calculate the sample points and weights, then map them, will use gaussxwab
#This function is more convenient because it does the mapping for us 

N1 = 8
x1, w1 = gaussxwab(N1, 0, x0)

#4.2 perform the integration, output calculated period value
s1 = 0.0
for i in range(N1):
    s1+= w1[i]*(f(x1[i], x0))

T1 = 4*s1

print('The calculated period for N=8 is', T1, 'seconds')

#5 Do the same for N=16

N2 = 16
x2, w2 = gaussxwab(N2, 0, x0)

s2 = 0.0
for i in range(N2):
    s2+= w2[i]*(f(x2[i], x0))

T2 = 4*s2

print('The calculated period for N=16 is', T2, 'seconds')


#6. Estimate error using Eqn 5.66 of the textbook, AND comparing with the classical value

error = s2 - s1 #Eqn. 5.66
print('The error estimate for these two N\'s is', error)

#Comparing with classical value

print('The fractional error for N=8 is:', (s1-(classical_T/4))/(classical_T/4))
print('The fractional error for N=16 is:', (s2-(classical_T/4))/(classical_T/4))



#Question 2b

print('Question 2b - see plots')

#Pseudocode for Q2b
#1. plot 4/gk for N=8 and N=16 at their sampling points
#2. plot 4wk/gk for N=8 and N=16 at their sampling points


#1. plot 4/gk for N=8 and N=16 at their sampling points

plt.figure()
plt.plot(x1,4/g(x1, x0), label='N=8')
plt.plot(x2,4/g(x2, x0), label='N=16')
plt.title('4/g(x) at the sampling points for N=8 and N=16')
plt.xlabel('Sampling point')
plt.ylabel('Integrand 4/gk')
plt.legend()
plt.show()


#2. plot 4wk/gk for N=8 and N=16 at their sampling points

plt.figure()
plt.plot(x1,4*w1/g(x1, x0), label='N=8')
plt.plot(x2,4*w2/g(x2, x0), label='N=16')
plt.title('4w/g(x) at the sampling points for N=8 and N=16')
plt.xlabel('Sampling point')
plt.ylabel('Weight value 4w/gk')
plt.legend()
plt.show()


#Question 2c
print('Question 2c - no code needed')


#Question 2d
print('Question 2d')

#Pseudocode for Q2d 
#1. Perform integration with gaussian quadrature for N= 200, copying code from previous parts
#2. Estimate error by finding fractional error and multiplying by 100%

#1. Perform integration with gaussian quadrature for N= 200

Nd = 200
xd, wd = gaussxwab(Nd, 0, x0)

sd = 0.0
for i in range(Nd):
    sd+= wd[i]*(f(xd[i], x0))

Td = 4*sd

print('The calculated period for N=200 is', Td, 'seconds')


#2. Estimate prercentage error by finding fractional error and multiplying by 100%

print('The percentage error for N=200 is:', ((sd-(classical_T/4))/(classical_T/4))*100, '%')


# Question 2e
print('Question 2e - see plots')

#Pseudocode for Q2e
#1. Initialize an array of x values for which we will perform the integration on
# using the value found in part c
#2. Initialize T_array, loop through x_array, calculate T at each point and save to T_array
#3. Plot



#1. Initialize an array of x values for which we will perform the integration on
# using the value found in part c
xc = c*np.sqrt(m/k) #metres

x_array = np.linspace(1, 8e8, 100) #the max value is 8e8 because 10*xc = 8.6e8


#2. Initialize T_array, loop through x_array, calculate T at each point and save to T_array

#Choose arbitrary Ne for guassian quadrature calculation
Ne = 200 #used this value because of its accuracy, and doesn't take to long to run each iteration

T_array = np.zeros(len(x_array))

for val in range(len(x_array)):
    xe, we = gaussxwab(Ne, 0, x_array[val])
    se = 0.0
    for i in range(Ne):
        se+= we[i]*(f(xe[i], x_array[val])) 
    T_array[val] = 4*se
    


#3. Plot T as a function of xc, also plot a zoomed in version to show classical limit    
plt.figure()
plt.plot(x_array, T_array)
plt.xlabel('x0')
plt.ylabel('Period')
plt.title('Period as a function of x0 in the range 1m < x < 10xc')
plt.show()


plt.figure()
plt.plot(x_array, T_array)
plt.xlabel('x0')
plt.ylabel('Period')
plt.title('Period as a function of x0 in the range 1m < x < 10xc, zoomed in')
plt.xlim(-0.2e8,0.8e8)
plt.show()
    
    
    
    


