#This code uses the shooting method and RK4 to find the bound states of hydrogen for a few different cases.
#Parts of this code (Question 2a) were taken from Newman's squarewell.py and modified. 

#Author: Milica Ivetic, October 2022

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

plt.rcParams['axes.titlesize'] = 13
plt.rcParams['axes.labelsize'] = 13
plt.rcParams['xtick.labelsize'] = 11
plt.rcParams['ytick.labelsize'] = 11

#Question 2a
print('Question 2a')

#Pseudocode for Question 2a
#1. Follow same format as Newman's squarewell.py
# i. Define constants and potential function
# ii. Define solve() function
# iii. Find energy using secant method


#Constants

m = 9.1094e-31 #mass of electron, kg
hbar = 1.0546e-34 #Planck's constant over 2*pi
e = 1.6022e-19 #Electron charge, Coulombs
a = 5.2918e-11 #Bohr radius, m
pi = np.pi
epsilon = 8.85418782e-12 #Vacuum permittivity, m-3 kg-1 s4 A2 units 

N = 1000 #Number of Runge-Kutta steps
h = 0.002*a #stepsize

#Potential function

def V(r):
    return -(e**2)/(4*pi*r*epsilon)

def f(w, r, E, n, l): #n is energy state
    R = (w[0])
    S = (w[1])
    
    #first-order ODEs calculated previously
    fR = S/(r**2)
    fS = R*(l*(l+1) + ((2*m*r**2)/(hbar**2))*(V(r) - E))
    return np.array([fR, fS], float)

#Calculate the wavefunction for a particular energy

def solve(E, n, l):
    
    R = 0.0
    S = 1.0
    w = np.array([R, S], float)
    R_array = np.array([R])
    
    
    #RK4 Method
    for i in np.arange(h, 20*a, h):
        k1 = h*f(w, i, E, n, l)
        k2 = h*f(w+0.5*k1, i+0.5*h, E, n, l)
        k3 = h*f(w+0.5*k2, i+0.5*h, E, n, l)
        k4 = h*f(w+k3, i+h, E, n, l)
        
        w += (k1 + 2*k2 + 2*k3 +k4)/6
        R_array = np.append(R_array, w[0])
        
        #print(R_array[-1], w[0]), ignore
        
    return w[0], R_array #R_array keeps track of R(r) at each timestep


#Main program to find the energy using the secant method

#Initializing eigen energies, depends on n
n= 1 #for now
l = 0 #for now
E1 = -15*e/(n**2)
E2 = -13*e/(n**2)

R2 = solve(E1, n, l)[0]

target = e/1000 #target energy convergence

while np.abs(E1 - E2) > target:
    R1, R2 = R2, solve(E2, n, l)[0]
    E1, E2 = E2, E2 - R2*(E2 - E1)/(R2 - R1)
    
    
print('E =', E2/e, 'eV. This was just a test run.')



#Question 2b
print('Question 2b')

#Pseudocode for question 2b
#1. Use code from previous part to calculate the energy for different values of n and l

#l = 0, n = 1

n= 1 #for now
l = 0 #for now
E1 = -15*e/(n**2)
E2 = -13*e/(n**2)

R2 = solve(E1, n, l)[0]

target = e/1000 #target energy convergence

while np.abs(E1 - E2) > target:
    R1, R2 = R2, solve(E2, n, l)[0]
    E1, E2 = E2, E2 - R2*(E2 - E1)/(R2 - R1)
    
    
print('E =', E2/e, 'eV for l = 0 and n = 1.')
print('The fractional error is', np.abs((E2/e - (-13.6))/(-13.6)))



#l = 0, n = 2

n= 2 #for now
l = 0 #for now
E1 = -15*e/(n**2)
E2 = -13*e/(n**2)

R2 = solve(E1, n, l)[0]


while np.abs(E1 - E2) > target:
    R1, R2 = R2, solve(E2, n, l)[0]
    E1, E2 = E2, E2 - R2*(E2 - E1)/(R2 - R1)
    
    
print('E =', E2/e, 'eV for l = 0 and n = 2.')
print('The fractional error is', np.abs((E2/e - (-3.4))/(-3.4)))



#l = 1, n = 2

n= 2 #for now
l = 1 #for now
E1 = -15*e/(n**2)
E2 = -13*e/(n**2)

R2 = solve(E1, n, l)[0]

while np.abs(E1 - E2) > target:
    R1, R2 = R2, solve(E2, n, l)[0]
    E1, E2 = E2, E2 - R2*(E2 - E1)/(R2 - R1)
    
    
print('E =', E2/e, 'eV for l = 1 and n = 2.')
print('The fractional error is', np.abs((E2/e - (-3.4))/(-3.4)))



#Question(2c)
print('Question 2c - see plots')

#Pseudocode for q2c
#1. Use scipy.integrate.cumtrapz and integrate over r for each case using trapezoidal method.
#2. Plot normalized eigenfunction for each case. 



#First case: l = 0, n = 1

n= 1 
l = 0 
E1 = -15*e/(n**2)
E2 = -13*e/(n**2)



R_sol = solve(E1, n, l)[1]
r_array = np.arange(h, 20*a, h)

norm = integrate.cumtrapz(np.abs(R_sol)**2)

plt.figure()
plt.yscale('log')
plt.plot(r_array/a, norm)
plt.title('Normalized eigenfunction R for l = 0, n = 1')
plt.xlabel('r/a [Bohr radius units]')
plt.ylabel('log(R)')
plt.show()




#Second case: l = 0, n = 2

n= 2 
l = 0 
E1 = -15*e/(n**2)
E2 = -13*e/(n**2)


R_sol = solve(E1, n, l)[1]
r_array = np.arange(h, 20*a, h)

norm = integrate.cumtrapz(np.abs(R_sol)**2)

plt.figure()
plt.yscale('log')
plt.plot(r_array/a, norm)
plt.title('Normalized eigenfunction R for l = 0, n = 2')
plt.xlabel('r/a [Bohr radius units]')
plt.ylabel('log(R)')
plt.show()


#Third case: l = 1, n = 2

n= 2 
l = 1 
E1 = -15*e/(n**2)
E2 = -13*e/(n**2)


R_sol = solve(E1, n, l)[1]
r_array = np.arange(h, 20*a, h)

norm = integrate.cumtrapz(np.abs(R_sol)**2)

plt.figure()
plt.yscale('log')
plt.plot(r_array/a, norm)
plt.title('Normalized eigenfunction R for l = 1, n = 2')
plt.xlabel('r/a [Bohr radius units]')
plt.ylabel('log(R)')
plt.show()



