#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from random import random, seed


# In[2]:


Tmax = 10.0
Tmin = 1e-8
tau = 2e4


# In[3]:


#define function 
def func(x,y):
    return x**2-np.cos(4*np.pi*x)+(y-1)**2


# In[4]:


#define initial conditions
ri = (2,2) #define inital coordinates r = (x,y)

x,y = ri


# In[5]:


# Function to generate two Gaussian random numbers from rutherford.py
sigma = 1.0
seed(10)
def gaussian():
    r = np.sqrt(-2*sigma*sigma*np.log(1-random()))
    theta = 2*np.pi*random()
    x = r*np.cos(theta)
    y = r*np.sin(theta)
    return x,y


# In[6]:


#Define arrays for solutions and add initial conditions
xvals_b = []
yvals_b = []

xvals_b.append(ri[0])
yvals_b.append(ri[1])


# In[7]:


#Simulaed annealing loop (adapted from salesman problem with suggestions from lab instructions)

t = 0
T = Tmax

#create time arrays for plotting
tvals_b = []
tvals_b.append(t)

while T>Tmin:
    # Cooling
    t += 1
    T = Tmax*np.exp(-t/tau)
    
    #Define random numbers drawn from a gaussain distribution
    dx, dy = gaussian()
    
    oldx = x
    oldy = y
    
    oldf = func(x,y)
    
    #apply "moves" x+dx and y+dy as noted in instructions 
    x += dx
    y += dy
    
    #calculate new f with shifted x,y values 
    f = func(x,y)
    
    #calculate change in function 
    deltaf = f-oldf
    
    # If the move is rejected, cancel the gaussian move 
    if random()>=np.exp(-deltaf/T):
        x = oldx
        y = oldy
        f = oldf
    
    #save points
    xvals_b.append(x)
    yvals_b.append(y)
    tvals_b.append(t)
    
print('Final value of (x,y) = ',(x,y))  


# In[9]:


##Plot values of (x, y) as a function of time
fig, ax = plt.subplots(1,2, dpi = 150)
fig.suptitle(r'$f(x,y)=x^2-cos(4\pi x)+(y-1)^2$')

ax[0].plot(tvals_b, xvals_b, '.')
ax[0].set_xlabel('t')
ax[0].set_ylabel('x position')

ax[1].plot(tvals_b, yvals_b, '.')
ax[1].set_xlabel('t')
ax[1].set_ylabel('y position')

ticks = [0,0,100000,200000,300000,400000]
ax[1].set_xticklabels(ticks, rotation=45, ha="right")
ax[0].set_xticklabels(ticks, rotation=45, ha="right")

#plt.xticks(rotation=45, ha="right")

fig.tight_layout()
fig.savefig('L11_Q1b.png', dpi = 100)
plt.show()


# In[15]:


#define function for part 1c
def func_c(x,y):
    if (0<x<50) and (-20<y<20):
        return np.cos(x)+np.cos(np.sqrt(2)*x)+np.cos(np.sqrt(3)*x)+(y-1)**2
    else: return 1E8


# In[16]:


#Define arrays for solutions and add initial conditions
xvals_c = []
yvals_c = []

xvals_c.append(ri[0])
yvals_c.append(ri[1])


# In[17]:


#Simulaed annealing loop (adapted from salesman problem with suggestions from lab instructions)

t = 0
T = Tmax

#create time arrays for plotting
tvals_c = []
tvals_c.append(t)

while T>Tmin:
    # Cooling
    t += 1
    T = Tmax*np.exp(-t/tau)
    
    #Define random numbers drawn from a gaussain distribution
    dx, dy = gaussian()
    
    oldx = x
    oldy = y
    
    oldfc = func_c(x,y)
    
    #apply "moves" x+dx and y+dy as noted in instructions
    x += dx
    y += dy
    
    #calculate new f with shifted x,y values 
    fc = func_c(x,y)
    
    #calculate changein function 
    deltafc = fc-oldfc
    
    # If the move is rejected, cancel the gaussian move 
    if random()>=np.exp(-deltafc/T):
        x = oldx
        y = oldy
        f = oldfc
    
    #save points
    xvals_c.append(x)
    yvals_c.append(y)
    tvals_c.append(t)
    
print('Final value of (x,y) = ',(x,y))  


# In[20]:


#Plot values of (x, y) as a function of time
fig, ax = plt.subplots(1,2, dpi = 150)
fig.suptitle(r'$f(x,y)= cos(x)+ cos(\sqrt{2}x)+cos(\sqrt{3}x)+(y-1)^2$')

ax[0].plot(tvals_c, xvals_c, '.')
ax[0].set_xlabel('t')
ax[0].set_ylabel('x position')

ax[1].plot(tvals_c, yvals_c, '.')
ax[1].set_xlabel('t')
ax[1].set_ylabel('y position')

ticks = [0,0,100000,200000,300000,400000]
ax[1].set_xticklabels(ticks, rotation=45, ha="right")
ax[0].set_xticklabels(ticks, rotation=45, ha="right")

fig.tight_layout()
fig.savefig('L11_Q1c.png', dpi = 100)
plt.show()


# In[ ]:




