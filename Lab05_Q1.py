# This code is an exercise in extracting a frequency spectrum from a time series 
#generated from a simulation of a relativistic particle on a spring.

#Author: Milica Ivetic, October 2022


import numpy as np
import matplotlib.pyplot as plt
from gaussxw import gaussxwab

plt.rcParams['axes.titlesize'] = 13
plt.rcParams['axes.labelsize'] = 13
plt.rcParams['xtick.labelsize'] = 11
plt.rcParams['ytick.labelsize'] = 11

#Question 1a
print('Question 1a  - see plots')


#Pseudocode for Question 1a

#1. Define constants, initialize time, position, and velocity arrays
#2. Using Euler-Cromer, simulate relativistic spring system from Lab3
#2a. euler-cromer for this system looks like:
     #v[i+1] = v[i] - dt*(k/m)*x[i]*(1 - (v[i]^2/c^2)^(3/2)
     #x[i+1] = x[i] + dt*v[i+1]
     
#3. Plot for each case


#1. Define constants, initialize time, position, and velocity arrays


#From lab 3, the relativistic spring system had the following properties

m = 1 #kg, mass
k = 12 #N/m, spring constant

c = 3.e8 #m/s, speed of light

#from Lab3, I know xc = c*sqrt(m/k)
xc = c*np.sqrt(m/k) #metres

#initial positions and velocities
x0_1 = 1 #metre
x0_2 = xc
x0_3 = 10*xc

v0 = 0 #m/s, starting from rest in all cases

#initializing arrays, need x and v for each case

dt = 0.001 #timestep
t = np.arange(0, 150, dt) #time array, in seconds

x1 = np.zeros(len(t))
v1 = np.zeros(len(t))

x2 = np.zeros(len(t))
v2 = np.zeros(len(t))

x3 = np.zeros(len(t))
v3 = np.zeros(len(t))

#applying initial conditions

x1[0] = x0_1
v1[0] = v0

x2[0] = x0_2
v2[0] = v0

x3[0] = x0_3
v3[0] = v0


#2. Using Euler-Cromer, simulate relativistic spring system from Lab3

#First case, x0 = 1m

for i in range(len(t) - 1):
    v1[i+1] = v1[i] - dt*(k/m)*x1[i]*(1 - v1[i]**2/c**2)**(3/2) 
    x1[i+1] = x1[i] + dt*v1[i+1]    



#Second case, x0 = xc

for i in range(len(t) - 1):
    v2[i+1] = v2[i] - dt*(k/m)*x2[i]*(1 - v2[i]**2/c**2)**(3/2) 
    x2[i+1] = x2[i] + dt*v2[i+1]  


#Third case, x0 = 10*xc

for i in range(len(t) - 1):
    v3[i+1] = v3[i] - dt*(k/m)*x3[i]*(1 - v3[i]**2/c**2)**(3/2) 
    x3[i+1] = x3[i] + dt*v3[i+1]  
    
    
#3. Plot for each case

plt.figure()
plt.title('Position vs. time of relativistic particle on spring for $x_0$ = 1m')
plt.plot(t, x1, c='blue')
plt.xlabel('Time [s]')
plt.ylabel('Position [m]')
plt.ylabel
plt.xlim(0,20)
plt.show()

plt.figure()
plt.title('Position vs. time of relativistic particle on spring for $x_0 = x_c$')
plt.plot(t, x2,c='orange')
plt.xlabel('Time [s]')
plt.ylabel('Position [m]')
plt.ylabel
plt.xlim(0,30)
plt.show()

plt.figure()
plt.title('Position vs. time of relativistic particle on spring for $x_0 = 10x_c$')
plt.plot(t, x3, c='green')
plt.xlabel('Time [s]')
plt.ylabel('Position [m]')
plt.ylabel
plt.show()


#Question 1b
print('Question 1b - see plot')

#Pseudocode for Question 1b
#1. For each case, calculate the Fourier transform using np.fft.fft
#2. Find the max values of each transformed array and calculate scaled arrays
#3. Plot the transforms on the same plot


#1. For each case, calculate the Fourier transform using np.fft.fft
x1_t = np.fft.rfft(x1) #Fourier transform
x2_t = np.fft.rfft(x2)
x3_t = np.fft.rfft(x3)


#2. Find the max values of each transformed array and calculate scaled arrays


x1t_max = np.amax(np.abs(x1_t))
x2t_max = np.amax(np.abs(x2_t))
x3t_max = np.amax(np.abs(x3_t))


x1t_scaled = np.abs(x1_t)/x1t_max
x2t_scaled = np.abs(x2_t)/x2t_max
x3t_scaled = np.abs(x3_t)/x3t_max



#Find x-axis (sampling frequencies) using np.fft.rfftfreq

sample_rate = 1/dt
x1_freq = np.fft.rfftfreq(x1.size, d=1/sample_rate)
x2_freq = np.fft.rfftfreq(x2.size, d=1/sample_rate)
x3_freq = np.fft.rfftfreq(x3.size,d=1/sample_rate)



#3. Plot the transforms on the same plot

plt.figure()
plt.title('Scaled Fourier transformed solutions for position vs. Frequency')
plt.plot(x1_freq, x1t_scaled, label = '$x_0$ = 1m')
plt.plot(x2_freq, x2t_scaled, label = '$x_0 = x_c$')
plt.plot(x3_freq, x3t_scaled, label = '$x_0 = 10x_c$')
plt.ylabel('Scaled Fourier transform')
plt.xlabel('Frequency [Hz]')
plt.xlim(0,1)
plt.legend()
plt.show()



#Question 1c
print('Question 1c - see plot')

#Pseudocode for Question 1c
#1. Calculate predicted period using Gaussian quadrature integration for each case
#1a. Define functions needed to calculate T from Eq.7 in Lab 3
#1b. Use loop from Newman, also used in Lab 3, repeat for each case
#2. Invert the calculated periods to get the frequencies for each case
#3. Plot vertical line with value 1/T on plot from before, for each case



#1. Calculate predicted period using Gaussian quadrature integration for each case
#1a. Define functions needed to calculate T from Eq.7 in Lab 3
def g(x, x0):
    diff = x0**2 - x**2
    m = 1 #kg, mass
    k = 12 #N/m, spring constant
    
    c = 3.e8 #m/s, speed of light    
    return c*((k*(diff)*(2*m*c**2+k*(diff)/2))/(2*(m*c**2+k*(diff)/2)**2))**(0.5)

def f(x, x0):
    return 1/g(x,x0)


#1b. Use loop from Newman, also used in Lab 3, calculate T, repeat for each case

#First case: x0 = 1m
N1 = 128
gauss_x1, gauss_w1 = gaussxwab(N1, 0, x0_1)

s1 = 0.0
for i in range(N1):
    s1+= gauss_w1[i]*(f(gauss_x1[i], x0_1))
T1 = 4*s1



#Second case: x0 = xc

gauss_x2, gauss_w2 = gaussxwab(N1, 0, x0_2)

s2 = 0.0
for i in range(N1):
    s2+= gauss_w2[i]*(f(gauss_x2[i], x0_2))
T2 = 4*s2



#Third case: x0 = 10xc

gauss_x3, gauss_w3 = gaussxwab(N1, 0, x0_3)

s3 = 0.0
for i in range(N1):
    s3+= gauss_w3[i]*(f(gauss_x3[i], x0_3))
T3 = 4*s3


#2. Invert the calculated periods to get the frequencies for each case

f1 = (1/T1)
f2 = (1/T2)
f3 = (1/T3)



#3. Plot vertical line with value 1/T on plot from before, for each case


print('frequencies',f1, f2, f3)

plt.figure()
plt.title('Scaled Fourier transformed solutions for position vs. Frequency')
plt.plot(x1_freq, x1t_scaled, label = '$x_0$ = 1m')
plt.plot(x2_freq, x2t_scaled, label = '$x_0 = x_c$')
plt.plot(x3_freq, x3t_scaled, label = '$x_0 = 10x_c$')
plt.ylabel('Scaled Fourier transform')
plt.xlabel('Frequency [Hz]')
plt.axvline(f1)
plt.axvline(f2)
plt.axvline(f3)
plt.xlim(0,1)
plt.legend()
plt.show()


