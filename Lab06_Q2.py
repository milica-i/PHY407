#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt


# In[2]:


#Q2(a) Pseudocode 
#i. define force from the negative gradient of the lj potential as a function of r(x,y) 
    #and return an 2d array of force components.
#ii. Define simulation and initial parameters for 16 partical system.
#iii. Define the acceleration from the F=ma as a function of x,y and return 2d 
    #array of acceleration components.
#iv. Define arrays for computed acceleration, velocity and position at each timestep 
    #in the given parameters dt=0.01, T=1000.
#v. Define initial acceleration and velocity from initial position conditions at t = 0.
#vi. Loop through time steps and apply the Verlet Theorem: find acceleration at time t
    #solve position at time t, plot position, update velocity.
#vii. Save position at each timestep in an array.


# In[3]:


def f(r, t):
    '''
    Defines x,y components of force from the negative gradient of lj potential 
    '''
    x,y = r
    #define f(r) components fx and fy from Q1a
    epsilon = 1
    m = 1
    sigma = 1
    
    fx = x*48*epsilon*((-1*sigma**6/(x**2+y**2)**7) + 2*(sigma**12/(x**2+y**2)**13))
    fy = y*48*epsilon*((-1*sigma**6/(x**2+y**2)**7) + 2*(sigma**12/(x**2+y**2)**13))
    return np.array([fx,fy]) 


# In[4]:


# Define simulation parameters
N         = 16    # Number of particles
tEnd      = 10.0   # time at which simulation ends
dt        = 0.01   # timestep
t = np.arange(0, 10, dt) #total time array

# Define initial paramters for 16 particles
Lx = 4.0
Ly = 4.0
dx = Lx/np.sqrt(N)
dy = Ly/np.sqrt(N)
x_grid = np.arange(dx/2, Lx, dx)
y_grid = np.arange(dy/2, Ly, dy)
xx_grid, yy_grid = np.meshgrid(x_grid, y_grid)

x_initial = xx_grid.flatten()
y_initial = yy_grid.flatten()

#create arrays for initial parameters 
posx = np.array(x_initial)
posy = np.array(y_initial)


# In[5]:


def accelerations(r, m, t):
    '''
    Def=ines x,y components of acceleration from force by newtons second law F=ma.
    '''
    #define acelleration arrays for each component x,y
    ax = np.zeros(posx.shape)
    ay = np.zeros(posy.shape)
    #determine spacing and acceleration between each of the 16 particles
    for i in range(len(x_initial)):
        for j in range(len(x_initial)):
            if i != j:
                dx = posx[i]-posx[j]
                dy = posy[i]-posy[j]
                
                a = f([dx,dy], t)
                ax[i] += a[0]/m
                ay[i] += a[1]/m
    return ax, ay


# In[6]:


accelerations((posx,posy), 1, 0)


# In[7]:


def position(posx,posy):
    '''
    Defines x,y components of position with the Vetlet thm and plots trajectories 
    of all 16 particles over 1000 iterations of timestep magnitude dt = 0.01.
    '''
    #define arrays for accelerations
    ax = np.zeros(posx.shape)
    ay = np.zeros(posy.shape)

    #define arrays for velocities
    vx = np.zeros(posx.shape)
    vy = np.zeros(posy.shape)

    m = 1
    
    #define initial acellerations
    a_i = accelerations((posx,posy),1,0)

    #define initial velocities
    vx = 0.5*dt * a_i[0] 
    vy = 0.5*dt * a_i[1]

    #define array for position data
    x = []
    y = []
    
    plt.figure(figsize = (10,8), dpi = 150)
    
    #Verlet Theorem loop 1000 steps 
    for t in np.arange(0,1000*dt, dt):
        #determine acelleration
        ax, ay = accelerations((posx,posy),1,t)
        
        #determine position 
        posx += dt*vx
        posy += dt*vy

        #plot positions
        plt.plot(posx,posy, 'k.')
        
        #add positions to position array
        x.append(posx)
        y.append(posy)
        
        #Update velocity with acceleration at time t 
        vx += dt*ax
        vy += dt*ay 
        
        
    plt.title('Trajectories of 16 Particles')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig('q2a.png', dpi = 150)
    plt.show()
    
    return x,y


# In[8]:


x,y = position(posx,posy)


# In[9]:


#Q2(b) Pseudocode 
#i. Define lj potential as a function of r(x,y) 
#ii. Define Energy as E = PE+KE as a function of r(x,y) 
#iii. Define arrays for computed energy at each timestep 
    #in the given parameters dt=0.01, T=1000.
#v. Define initial acceleration and velocity from initial position conditions at t = 0.
#vi. Loop through time steps and apply the Verlet Theorem: find acceleration at time t
    #solve position at time t, determine energy, plot energy, update velocity.
#vii. Save energy at each timestep in an array.


# In[10]:


def lj_potential(r,t):
    '''
    Defines lj potential.
    '''
    x,y = r
    epsilon = 1
    sigma = 1
    return 4*epsilon*((sigma/(x**2+y**2))**12 - (sigma/(x**2+y**2))**6)

def energy_comp(r,v,t):
    '''
    Defines total system energy - sum of kinetic and potential energy.
    '''
    m = 1
    KE = 0.5*m*v**2
    PE = lj_potential(r,t)
    E = KE+PE
    return E

def energy(posx,posy):
    '''
    Defines energy at each point in the time array with the Vetlet thm and plots energy 
    of all 16 particles over 1000 iterations of timestep magnitude dt = 0.01.
    '''
    #define arrays for energy
    E = np.zeros(posx.shape)
    
    #define initial acellerations
    a_i = accelerations((posx,posy),1,0)
    
    #define initial velocities
    vx = 0.5*dt * a_i[0] 
    vy = 0.5*dt * a_i[1]
    
    #define energies
    E_i = energy_comp((posx,posy), np.sqrt(vx**2+vy**2),0)

    #Verlet Theorem loop 1000 steps 
    for t in np.arange(0,1000*dt, dt):
        #determine acelleration
        ax, ay = accelerations((posx,posy),1,t)
        
        #determine position 
        posx += dt*vx
        posy += dt*vy
        
        v = np.sqrt(vx**2+vy**2)
        
        E +=  energy_comp((posx,posy), v,t)
        for n in range(16):
            plt.plot(t,E[n], '.')
        #Update velocity with acceleration at time t 
        vx += dt*ax
        vy += dt*ay 
    plt.show()
    return E
    


# In[ ]:


E = energy(posx,posy)


# In[ ]:




