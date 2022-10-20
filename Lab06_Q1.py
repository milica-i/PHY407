#This code considers a 2D molecular dynamics simulation of N = 2 molecules under 
#the influence of only the Lennard-Jones potential.
#Positions are updated using Verlet method.

#Author: Milica Ivetic, October 2022


import numpy as np
import matplotlib.pyplot as plt


#Question 1b
print('Question 1b')

#Pseudocode for Question 1b:
#1. Initialize time array, velocities, x and y arrays, and function f we are working with
#2. Implement Verlet method on both particles for 3 different sets of initial conditions
#2a: Verlet method is as follows:
#For the first step: v[t+dt/2] = v[t] + (dt/2)*f(r[t], t)
# The following steps: r[t + dt] = r[t] + dt*v[t+dt/2]
                              #k = dt*f(r[t+dt], t+dt)
                              #v[t+h] = v[t+dt/2] + (1/2)*k
                              #v[t+3*dt/2] = v[t + dt/2] + k
# Where r is x or y depending on which component we are calculating
#3. Plot trajectories for each set of initial conditions


#1. Initialize time array, velocities, x and y arrays, and functions fx and fy we are working with

def fx(x1, x2, y1, y2):
                              w = x2 - x1
                              z = y2 - y1
                              r = np.sqrt(w**2 + z**2)
                              
                              return -4*(12/(r)**14 - 6/(r)**8)*w


def fy(x1, x2, y1, y2):
                              w = x2 - x1
                              z = y2 - y1
                              r = np.sqrt(w**2 + z**2)
                              
                              return -4*(12/(r)**14 - 6/(r)**8)*z


v0 = 0 #m/s, the same initial velocities for all cases

dt = 0.01 #s
t = np.arange(0, 1, dt)


#First set of initial conditions, called set "a"
x1a0 = 4.
y1a0 = 4.

x2a0 = 5.2
y2a0 = 4.


#Second set of initial conditions, called set "b"
x1b0 = 4.5
y1b0 = 4.

x2b0 = 5.2
y2b0 = 4.


#Third set of initial conditions, called set "c"
x1c0 = 2.
y1c0 = 3.

x2c0 = 3.5
y2c0 = 4.4

#x, y, vx, and vy arrays for each set of ICs

#a
x1a = np.array([], dtype=object)
y1a = np.array([], dtype=object)
x2a = np.array([], dtype=object)
y2a = np.array([], dtype=object)

vx1a = np.array([], dtype=object)
vy1a = np.array([], dtype=object)
vx2a = np.array([], dtype=object)
vy2a = np.array([], dtype=object)


#b
x1b = np.array([], dtype=object)
y1b = np.array([], dtype=object)
x2b = np.array([], dtype=object)
y2b = np.array([], dtype=object)

vx1b = np.array([], dtype=object)
vy1b = np.array([], dtype=object)
vx2b = np.array([], dtype=object)
vy2b = np.array([], dtype=object)


#c
x1c = np.array([], dtype=object)
y1c = np.array([], dtype=object)
x2c = np.array([], dtype=object)
y2c = np.array([], dtype=object)

vx1c = np.array([], dtype=object)
vy1c = np.array([], dtype=object)
vx2c = np.array([], dtype=object)
vy2c = np.array([], dtype=object)


#initialize

#a
#x1a.append(x1a0)
x1a = np.append(x1a, x1a0)
print('dfd',x1a)
y1a = np.append(y1a,y1a0)
x2a = .append(x2a0)
y2a.append(y2a0)

vx1a.append(v0)
vy1a.append(v0)
vx2a.append(v0)
vy2a.append(v0)


#b
x1b.append(x1b0)
y1b.append(y1b0)
x2b.append(x2b0)
y2b.append(y2b0)

vx1b.append(v0)
vy1b.append(v0)
vx2b.append(v0)
vy2b.append(v0)

#c
x1c.append(x1c0)
y1c.append(y1c0)
x2c.append(x2c0)
y2c.append(y2c0)

vx1c.append(v0)
vy1c.append(v0)
vx2c.append(v0)
vy2c.append(v0)


#2. Implement Verlet method on both particles for 3 different sets of initial conditions

#Set a

kx1a = []
ky1a = []
kx2a = []
ky2a = []

#First time step:

for i in range(len(t)-1):
                              if i == 0:
                                                            print('dsfds',vx1a[-1], vx2a[-1], vy1a[-1], vy2a[-1])
                                                            #v[t+dt/2] = v[t] + (dt/2)*f(r[t], t)
                                                            vx1a.append(vx1a + (dt/2)*fx(x1a[-1], x2a[-1], y1a[-1], y2a[-1]))
                                                            vy1a.append(vy1a + (dt/2)*fy(x1a[-1], x2a[-1], y1a[-1], y2a[-1]))
                                                            vx2a.append(vx1a - (dt/2)*fx(x1a[-1], x2a[-1], y1a[-1], y2a[-1]))
                                                            vy2a.append(vy1a - (dt/2)*fy(x1a[-1], x2a[-1], y1a[-1], y2a[-1]))
                                                            print('func',fy(x1a[-1], x2a[-1], y1a[-1], y2a[-1]))
                                                            print('func',fx(x1a[-1], x2a[-1], y1a[-1], y2a[-1]))
                                                            print('dsfds',vx1a[-1], vx2a[-1], vy1a[-1], vy2a[-1])
                              else:
                                                            print(i)
                                                            #r[t + dt] = r[t] + dt*v[t+dt/2]
                                                            #k = dt*f(r[t+dt], t+dt)
                                                            #v[t+h] = v[t+dt/2] + (1/2)*k -> NOT NEEDED FOR TIMESTEPPING
                                                            #v[t+3*dt/2] = v[t + dt/2] + k   
                                                            
                                                            x1a.append(x1a[-1] + dt*vx1a[-1])
                                                            y1a.append(y1a[-1] + dt*vy1a[-1])
                                                            x2a.append(x2a[-1] + dt*vx2a[-1])
                                                            y2a.append(y2a[-1] + dt*vy2a[-1])
                                                            
                                                            print(x1a[-1], x2a[-1], y1a[-1], y2a[-1])
                                                            print('cvxcxv',fx(x1a[-1], x2a[-1], y1a[-1], y2a[-1]))
                                                            kx1a.append(dt*fx(x1a[-1], x2a[-1], y1a[-1], y2a[-1]))
                                                            ky1a.append(dt*fy(x1a[-1], x2a[-1], y1a[-1], y2a[-1]))
                                                            kx2a.append(dt*(-1)*fx(x1a[-1], x2a[-1], y1a[-1], y2a[-1]))
                                                            ky2a.append(dt*(-1)*fy(x1a[-1], x2a[-1], y1a[-1], y2a[-1]))
                                                            
                                                            vx1a.append(vx1a[-1] + kx1a[-1])
                                                            vy1a.append(vy1a[-1] + ky1a[-1])
                                                            vx2a.append(vx2a[-1] + kx2a[-1])
                                                            vy2a.append(vy2a[-1] + ky2a[-1])
                              
                              
print(x1a, y1a)
                                                            
                                                            
                                                       
                                                            

