#!/usr/bin/env python
# coding: utf-8

# In[1]:


'''
Written by Madeline Nardin Oct. 2022
'''
import numpy as np 
import matplotlib.pyplot as plt
import time 


# In[2]:


#Q1a Pseudocode 
#i. Use  Nico Grisouard Solution to Newman 8.8 to determine fixed timestep orbital 
    #dynamics.
#ii. Time fixed method execution duration for part b
#Adaptive method
#iii. Define error conditions and empty lists for determined x,vx,y,vy, t, and h at 
    #each point in time.
#iv. Define function for rho from textbook definition
#v. Define inital time value and start timer
#vi. Preform RK4 method with adaptive step size: first step with inital r and step size
    #h, second step with r computed in first step and step size h, third step with initial
    #r and step size 2h.
#vii. Compute rho and impose rho>1 condition to ensure accuracy dominates target accuracy
#viii. If accuracy condition is met save components of r to the respective lists and step time, 
    #if not repeat with next value t.
#ix. Stop timer and print execution duration.
#x. Plot both adaptive and fixed time step method on the same plot.


# In[3]:


# Note the following code has been adapted from Solution to Newman 8.8, 
#Space garbage.Author: Nico Grisouard, Univ. of Toronto --------------|
def rhs(r,t):
    """ The right-hand-side of the equations
    INPUT:
    r = [x, vx, y, vy] are floats (not arrays)
    note: no explicit dependence on time
    OUTPUT:
    1x2 numpy array, rhs[0] is for x, rhs[1] is for vx, etc"""
    M = 10.
    L = 2.

    x = r[0]
    vx = r[1]
    y = r[2]
    vy = r[3]

    r2 = x**2 + y**2
    Fx, Fy = - M * np.array([x, y], float) / (r2 * np.sqrt(r2 + .25*L**(2)/4))
    return np.array([vx, Fx, vy, Fy], float)

#Define iteration conditions
ti = 0.0
tf = 10.0
N = 10000 
h = (tf-ti)/N

#Define time array 
tpoints = np.arange(ti, tf, h)

#define arrays for determined x,vx,y,vy at each point in time 
xpoints = []
vxpoints = []  # the future dx/dt
ypoints = []
vypoints = []  # the future dy/dt
    
#Define initial conditions for x,vx,y,vy
r = np.array([1., 0., 0., 1.], float)

# Preform 4th order Runge-Kutta (RK4)
#Start timer
fixed_i = time.time()
for t in tpoints:
    #set initial conditions
    xpoints.append(r[0])
    vxpoints.append(r[1])
    ypoints.append(r[2])
    vypoints.append(r[3])
    
    #note for non-adaptive (fixed) time step rhs is independant of time 
    k1 = h*rhs(r,0)  
    k2 = h*rhs(r + 0.5*k1,0)  
    k3 = h*rhs(r + 0.5*k2,0)
    k4 = h*rhs(r + k3,0)
    r += (k1 + 2*k2 + 2*k3 + k4)/6
#Stop timer
fixed_f = time.time()
#Print execution duration
print('Fixed time step computation duration = ', fixed_f-fixed_i, 's')
#----------------------------------------------------------------------|


# In[5]:


# Set initial parameters for adaptive step size method 
delta = 10**(-6)

#define arrays for determined x,vx,y,vy and time step size (h) at each point in time 
adaptive_tpoints = []
adaptive_hpoints = []
adaptive_xpoints = []
adaptive_vxpoints = []
adaptive_ypoints = []
adaptive_vypoints = []

#Define accuracy ratio rho 
def rho(h, delta, x1, x2, y1, y2):
    ''' rho is the ratio of the target accuracy (product of target accuracy and step size) 
    and the actual accuracy, given by rho = delta*h/sqrt(epsx**2+epsy**2) such that epsx 
    and epsy are the error terms in the x and y component given by epsi = (i1-i2)/30.
    ''' 
    epsx = (x1 - x2)/30
    epsy = (y1 - y2)/30

    d = np.sqrt(epsx**2 + epsy**2)
    return h*delta/d

#Define initial time
t = ti

#Start timer
adaptive_i = time.time()

# Preform RK4 method with adaptive step size
while t < tf:
    #First step of size h starting from inital conditions
    k1 = h * rhs(r, t)
    k2 = h * rhs(r + 0.5 * k1, t + 0.5 * h)
    k3 = h * rhs(r + 0.5 * k2, t + 0.5 * h)
    k4 = h * rhs(r + k3, t + h)
    r0 = r + (k1 + 2 * k2 + 2 * k3 + k4) / 6

    #Second step of size h starting from r0 computed in first step with single step size increase of h
    tstepped = t + h
    k1 = h * rhs(r0, tstepped)
    k2 = h * rhs(r0 + 0.5 * k1, tstepped + 0.5 * h)
    k3 = h * rhs(r0 + 0.5 * k2, tstepped + 0.5 * h)
    k4 = h * rhs(r0 + k3, tstepped + h)
    r1 = r0 + (k1 + 2 * k2 + 2 * k3 + k4) / 6

    #Third step of size 2h starting from initial conditions with step size increase of 2h
    hstepped = 2 * h
    k1 = hstepped * rhs(r, t)
    k2 = hstepped * rhs(r + 0.5 * k1, t + 0.5 * hstepped)
    k3 = hstepped * rhs(r + 0.5 * k2, t + 0.5 * h)
    k4 = hstepped * rhs(r + k3, t + hstepped)
    r2 = r + (k1 + 2 * k2 + 2 * k3 + k4) / 6

    # Compute rho
    rho_val = rho(h, delta, r1[0], r2[0], r1[2], r2[2])

    # Compute new h
    h *= rho_val**(0.25)
    
    #If rho>1 actual accuracy dominates thus we save values, if not the process repeats for the next
    #value of t.
    if rho_val > 1.0:
        # Append time and stepsize to the respective lists 
        adaptive_tpoints.append(t)
        adaptive_hpoints.append(h)

        #Step time
        t += 2 * h

        #Save r as a vector 
        r = r1

        #Save components of r to the respective lists
        adaptive_xpoints.append(r[0])
        adaptive_vxpoints.append(r[1])
        adaptive_ypoints.append(r[2])
        adaptive_vypoints.append(r[3])
        
#Stop timer
adaptive_f= time.time()
#Print execution duration
print('Adaptive time step computation duration = ', adaptive_f-adaptive_i, 's')


# In[6]:


#Plot both adaptive and fixed time step method on the same plot
plt.figure(figsize = (10,8), dpi = 150)
plt.plot(xpoints, ypoints, 'k-', label = 'Fixed Step')
plt.plot(adaptive_xpoints,adaptive_ypoints, 'r.',label = 'Adapted Step')
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.title('Trajectory of a ball bearing around a space rod.')
plt.axis('equal')
plt.grid()
plt.tight_layout()
plt.savefig('Lab071a.png', dpi=150)
plt.legend()


# In[7]:


#Plot time step size as a funtion of time from saved values in adaptive step method
plt.figure(dpi = 150)
plt.plot(adaptive_tpoints, adaptive_hpoints, 'k.-')
plt.xlabel("$t$")
plt.ylabel("Time Step Size")
plt.title('Time Step Size as a Function of Time')
plt.ylim(0,0.07)
plt.grid()
plt.tight_layout()
plt.savefig('Lab071c.png', dpi=150)


# In[ ]:




