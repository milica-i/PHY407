# The following code explores roundoff error and its effects on a polynomial.
# Author: Milica Ivetic

import numpy as np
import matplotlib.pyplot as plt

#The functions we are considering in this question (p & q):

def p(u):
    return (1-u)**8

def q(u):
    return 1 - 8*u + 28*(u**2) - 56*(u**3) + 70*(u**4) - 56*(u**5) + 28*(u**6) - 8*(u**7) + u**8

# Question 4a

# Plotting p and q on the same graph for values of u that are very close to 1

u_array = np.linspace(0.98, 1.02, 500)

plt.figure()
plt.plot(u_array, q(u_array), label = 'q(u)')
plt.plot(u_array, p(u_array), label = 'p(u)')
plt.xlabel('u')
plt.ylabel('Function outputs')
plt.title('p(u) and q(u) for values of u very close to 1')
plt.legend()
plt.show()



# Question 4b

print('Question 4b')

#plotting p - q
plt.figure()
plt.plot(u_array, p(u_array) - q(u_array))
plt.title('p(u) - q(u) for values of u close to 1')
plt.xlabel('u')
plt.ylabel('p(u) - q(u)')
plt.show()

#plotting the histogram of p - q
plt.figure()
plt.hist(p(u_array) - q(u_array), bins=20)
plt.axvline(x= np.mean(p(u_array) - q(u_array)), c='red', label='Mean')
plt.title('Histogram of p(u) - q(u)')
plt.ylabel('Frequency')
plt.xlabel('p(u) - q(u)')
plt.legend()
plt.show()

# Calculating standard deviation of this distribution using numpy

std = np.std(p(u_array) - q(u_array))
print('Using np.std, the standard deviation is', std)


#Calculating the estimate using equation 3 of the lab

C = 1 # as per the announcement, in our case C=1
N = len(u_array) #Number of terms, which is the length of the array of u values

# root mean squared calculation, take the mean of p-q squared, and square root that
rms = np.sqrt(np.mean((p(u_array) - q(u_array))**2))

sigma = C*np.sqrt(N)*rms 

print('Using Equation 3 of the lab, sigma is', sigma)





# Question 4c

print('Question 4c')

#Using equation 4 from the lab
# For fraction error to be 100%, the LHS should be 1, we want to show
# that the RHS is also =~1 for u=~0.980

#Defining a new set of u values to consider
u_array_2 = np.linspace(0.982, 0.990, 500)
N2 = len(u_array_2)

#new root mean squared value

rms_2 = np.sqrt(np.mean((p(u_array_2) - q(u_array_2))**2))

#new mean value
mean = np.mean(p(u_array_2) - q(u_array_2))

#Calculating fractional error
frac_error = (C*rms_2)/(np.sqrt(N2)*mean)
print('Using 0.982 < u < 0.990, the fractional error is', frac_error)




# Verify using abs(p-q)/abs(p)
#Defining a new set of u values to consider
u_array_3 = np.linspace(0.980, 0.984, 1000)

verification = np.abs(p(u_array_3) - q(u_array_3))/np.abs(p(u_array_3))

#plot this quantity
plt.figure()
plt.title('abs(p-q)/abs(p) for 0.980 < u < 0.984')
plt.xlabel('u')
plt.ylabel('abs(p-q)/abs(p)')
plt.plot(u_array_3, verification)
plt.show()



# Question 4d
print('Question 4d')

#Back to using the first u_array

#define new function f

def f(u):
    return (u**8)/((u**4)*(u**4))

std_f = np.std(f(u_array))

print('Standard deviation using numpy gives', std_f)

plt.figure()
plt.title('f(u)-1 for values of u close to 1')
plt.ylabel('f(u) - 1')
plt.xlabel('u')
plt.plot(u_array, f(u_array) - 1)
plt.ylim(-1e-15,1e-15)
plt.show()


#Comparing calculated standard deviation with equation 4.5 of the textbook

C2 = 1e-16
error_d = np.sqrt(2)*C2*f(u_array)

print('Error estimate gives', error_d[-1])

