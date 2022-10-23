#This code considers a 2D molecular dynamics simulation of N = 2 molecules under 
#the influence of only the Lennard-Jones potential.
#Positions are updated using Verlet method.

#Author: Milica Ivetic, October 2022


import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['axes.titlesize'] = 13
plt.rcParams['axes.labelsize'] = 13
plt.rcParams['xtick.labelsize'] = 11
plt.rcParams['ytick.labelsize'] = 11


#Question 1b
print('Question 1b')

#Pseudocode for Question 1b:
#1. Initialize time array, velocities, x and y arrays, and functions we are working with
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

#x-component
def fx(x1, x2, y1, y2):
                              w = x2 - x1
                              z = y2 - y1
                              r = np.sqrt(w**2 + z**2)
                              
                              return -4*(12/(r)**14 - 6/(r)**8)*w

#y-component
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
x1a = np.append(x1a, x1a0)
y1a = np.append(y1a, y1a0)
x2a = np.append(x2a, x2a0)
y2a = np.append(y2a, y2a0)

vx1a = np.append(vx1a, v0)
vy1a = np.append(vy1a, v0)
vx2a = np.append(vx2a, v0)
vy2a = np.append(vy2a, v0)


#b
x1b = np.append(x1b, x1b0)
y1b = np.append(y1b, y1b0)
x2b = np.append(x2b, x2b0)
y2b = np.append(y2b, y2b0)

vx1b = np.append(vx1b, v0)
vy1b = np.append(vy1b, v0)
vx2b = np.append(vx2b, v0)
vy2b = np.append(vy2b, v0)

#c
x1c = np.append(x1c, x1c0)
y1c = np.append(y1c, y1c0)
x2c = np.append(x2c, x2c0)
y2c = np.append(y2c, y2c0)

vx1c = np.append(vx1c, v0)
vy1c = np.append(vy1c, v0)
vx2c = np.append(vx2c, v0)
vy2c = np.append(vy2c, v0)


#2. Implement Verlet method on both particles for 3 different sets of initial conditions

#Set a

kx1a = np.array([], dtype=object)
ky1a = np.array([], dtype=object)
kx2a = np.array([], dtype=object)
ky2a = np.array([], dtype=object)

#First time step:

for i in range(len(t)):
                              if i == 0:
                                                            #v[t+dt/2] = v[t] + (dt/2)*f(r[t], t)
                                                            vx1a = np.append(vx1a, vx1a + (dt/2)*fx(x1a[-1], x2a[-1], y1a[-1], y2a[-1]))
                                                            vy1a = np.append(vy1a, vy1a + (dt/2)*fy(x1a[-1], x2a[-1], y1a[-1], y2a[-1]))
                                                            vx2a = np.append(vx2a, vx2a - (dt/2)*fx(x1a[-1], x2a[-1], y1a[-1], y2a[-1]))
                                                            vy2a = np.append(vy2a, vy2a - (dt/2)*fy(x1a[-1], x2a[-1], y1a[-1], y2a[-1]))
                              
                              else:
                                                            #r[t + dt] = r[t] + dt*v[t+dt/2]
                                                            #k = dt*f(r[t+dt], t+dt)
                                                            #v[t+h] = v[t+dt/2] + (1/2)*k -> NOT NEEDED FOR TIMESTEPPING
                                                            #v[t+3*dt/2] = v[t + dt/2] + k   
                                                            
                                                            x1a = np.append(x1a, x1a[-1] + dt*vx1a[-1])
                                                            y1a = np.append(y1a, y1a[-1] + dt*vy1a[-1])
                                                            x2a = np.append(x2a, x2a[-1] + dt*vx2a[-1])
                                                            y2a = np.append(y2a, y2a[-1] + dt*vy2a[-1])
                                                            
                                                            kx1a = np.append(kx1a, dt*fx(x1a[-1], x2a[-1], y1a[-1], y2a[-1]))
                                                            ky1a = np.append(ky1a, dt*fy(x1a[-1], x2a[-1], y1a[-1], y2a[-1]))
                                                            kx2a = np.append(kx2a, dt*(-1)*fx(x1a[-1], x2a[-1], y1a[-1], y2a[-1]))
                                                            ky2a = np.append(ky2a, dt*(-1)*fy(x1a[-1], x2a[-1], y1a[-1], y2a[-1]))
                                                            
                                                            vx1a = np.append(vx1a, vx1a[-1] + kx1a[-1])
                                                            vy1a = np.append(vy1a, vy1a[-1] + ky1a[-1])
                                                            vx2a = np.append(vx2a, vx2a[-1] + kx2a[-1])
                                                            vy2a = np.append(vy2a, vy2a[-1] + ky2a[-1])
                              
                                                                                                                 
#Set b

kx1b = np.array([], dtype=object)
ky1b = np.array([], dtype=object)
kx2b = np.array([], dtype=object)
ky2b = np.array([], dtype=object)

#First time step:

for i in range(len(t)):
                              if i == 0:
                                                          
                                                            #v[t+dt/2] = v[t] + (dt/2)*f(r[t], t)
                                                            vx1b = np.append(vx1b, vx1b + (dt/2)*fx(x1b[-1], x2b[-1], y1b[-1], y2b[-1]))
                                                            vy1b = np.append(vy1b, vy1b + (dt/2)*fy(x1b[-1], x2b[-1], y1b[-1], y2b[-1]))
                                                            vx2b = np.append(vx2b, vx2b - (dt/2)*fx(x1b[-1], x2b[-1], y1b[-1], y2b[-1]))
                                                            vy2b = np.append(vy2b, vy2b - (dt/2)*fy(x1b[-1], x2b[-1], y1b[-1], y2b[-1]))
                              
                              else:
                                                            #r[t + dt] = r[t] + dt*v[t+dt/2]
                                                            #k = dt*f(r[t+dt], t+dt)
                                                            #v[t+h] = v[t+dt/2] + (1/2)*k -> NOT NEEDED FOR TIMESTEPPING
                                                            #v[t+3*dt/2] = v[t + dt/2] + k   
                                                            
                                                            x1b = np.append(x1b, x1b[-1] + dt*vx1b[-1])
                                                            y1b = np.append(y1b, y1b[-1] + dt*vy1b[-1])
                                                            x2b = np.append(x2b, x2b[-1] + dt*vx2b[-1])
                                                            y2b = np.append(y2b, y2b[-1] + dt*vy2b[-1])
                                                            
                                                            
                                                            kx1b = np.append(kx1b, dt*fx(x1b[-1], x2b[-1], y1b[-1], y2b[-1]))
                                                            ky1b = np.append(ky1b, dt*fy(x1b[-1], x2b[-1], y1b[-1], y2b[-1]))
                                                            kx2b = np.append(kx2b, dt*(-1)*fx(x1b[-1], x2b[-1], y1b[-1], y2b[-1]))
                                                            ky2b = np.append(ky2b, dt*(-1)*fy(x1b[-1], x2b[-1], y1b[-1], y2b[-1]))
                                                            
                                                            vx1b = np.append(vx1b, vx1b[-1] + kx1b[-1])
                                                            vy1b = np.append(vy1b, vy1b[-1] + ky1b[-1])
                                                            vx2b = np.append(vx2b, vx2b[-1] + kx2b[-1])
                                                            vy2b = np.append(vy2b, vy2b[-1] + ky2b[-1])
                                                            
                                                            

#Set c
                                                            
kx1c = np.array([], dtype=object)
ky1c = np.array([], dtype=object)
kx2c = np.array([], dtype=object)
ky2c = np.array([], dtype=object)
                                                            
#First time step:
                                                            
for i in range(len(t)):
                              if i == 0:
                                                                                                                      
                                                            #v[t+dt/2] = v[t] + (dt/2)*f(r[t], t)
                                                            vx1c = np.append(vx1c, vx1c + (dt/2)*fx(x1c[-1], x2c[-1], y1c[-1], y2c[-1]))
                                                            vy1c = np.append(vy1c, vy1c + (dt/2)*fy(x1c[-1], x2c[-1], y1c[-1], y2c[-1]))
                                                            vx2c = np.append(vx2c, vx2c - (dt/2)*fx(x1c[-1], x2c[-1], y1c[-1], y2c[-1]))
                                                            vy2c = np.append(vy2c, vy2c - (dt/2)*fy(x1c[-1], x2c[-1], y1c[-1], y2c[-1]))
                                                                                          
                              else:
                                                            #r[t + dt] = r[t] + dt*v[t+dt/2]
                                                            #k = dt*f(r[t+dt], t+dt)
                                                            #v[t+h] = v[t+dt/2] + (1/2)*k -> NOT NEEDED FOR TIMESTEPPING
                                                            #v[t+3*dt/2] = v[t + dt/2] + k   
                                                                                                                        
                                                            x1c = np.append(x1c, x1c[-1] + dt*vx1c[-1])
                                                            y1c = np.append(y1c, y1c[-1] + dt*vy1c[-1])
                                                            x2c = np.append(x2c, x2c[-1] + dt*vx2c[-1])
                                                            y2c = np.append(y2c, y2c[-1] + dt*vy2c[-1])
                                                                                                                                                                     
                                                            kx1c = np.append(kx1c, dt*fx(x1c[-1], x2c[-1], y1c[-1], y2c[-1]))
                                                            ky1c = np.append(ky1c, dt*fy(x1c[-1], x2c[-1], y1c[-1], y2c[-1]))
                                                            kx2c = np.append(kx2c, dt*(-1)*fx(x1c[-1], x2c[-1], y1c[-1], y2c[-1]))
                                                            ky2c = np.append(ky2c, dt*(-1)*fy(x1c[-1], x2c[-1], y1c[-1], y2c[-1]))
                                                                                                                        
                                                            vx1c = np.append(vx1c, vx1c[-1] + kx1c[-1])
                                                            vy1c = np.append(vy1c, vy1c[-1] + ky1c[-1])
                                                            vx2c = np.append(vx2c, vx2c[-1] + kx2c[-1])
                                                            vy2c = np.append(vy2c, vy2c[-1] + ky2c[-1])

          
                                           
#3. Plot trajectories for each set of initial conditions


#set a


plt.figure()
plt.plot(x1a, y1a, '.', label = '1')
plt.plot(x2a, y2a, '.', label = '2')
plt.title('Trajectories of both particles: first set of initial conditions')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.show()


#x vs. time of both particles, don't need to submit
#OSCILLATION IN X VS. T FOR BOTH PARTICLES
#plt.figure()
#plt.scatter(t, x1a)
#plt.scatter(t, x2a)
#plt.show()

#plt.figure()
#plt.scatter(t, y1a)
#plt.scatter(t, y2a)
#plt.show()


#set b

plt.figure()
plt.plot(x1b, y1b, '.', label = '1')
plt.plot(x2b, y2b, '.', label = '2')
plt.title('Trajectories of both particles: second set of initial conditions')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.show()

#x vs. time of both particles, don't need to submit
#Particles fly away from each other
#plt.figure()
#plt.title('set b')
#plt.scatter(t, x1b)
#plt.scatter(t, x2b)
#plt.ylabel('x')
#plt.show()

#plt.figure()
#plt.title('set b')
#plt.scatter(t, y1b)
#plt.scatter(t, y2b)
#plt.ylabel('y')
#plt.show()



#set c

plt.figure()
plt.plot(x1c, y1c, '.', label = '1')
plt.plot(x2c, y2c, '.', label = '2')
plt.title('Trajectories of both particles: third set of initial conditions')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.show()

#x vs. time of both particles, don't need to submit
#Particles start attracting 
#plt.figure()
#plt.title('set c')
#plt.scatter(t, x1c)
#plt.scatter(t, x2c)
#plt.show()

#plt.figure()
#plt.scatter(t, y1c)
#plt.scatter(t, y2c)
#plt.show()

