#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt


# In[2]:


#Define Constants 
g = 9.81 #[m/s^2]
A = 0.002
mu = 0.5 #[m]
sigma = 0.05 #[m]
nb = 0.0
H = 0.01 #[m]

#Define position bounds 
xi = 0.0 #[m]
xf = 1.0 #[m]

J = 50 #number of divisions in grid
dx = 0.02 #position step size [m]
x = np.arange(xi,xf,dx)

#Define time bounds
dt = 0.01 #time step size [s]
ti = 0.0 #[s]
tf = 4.0 #[s]

#Define times of interest 
t1 = 0.0 #[s]
t2 = 1.0 #[s]
t3 = 4.0 #[s]

eps = dt*1e-3


# In[3]:


#Define flux conserved form vector
def flux_cons(u,n,H):
    Fu = 0.5*(u**2) + g*n
    Fn = (n-nb)*u
    return Fu, Fn


# In[4]:


#Define initial conditions 
u0 = np.zeros(len(x))
eta0 = H+A*np.exp(-(x-mu)**2/(sigma)**2) - np.mean(A*np.exp(-(x-mu)**2/(sigma)**2))

#Define boundary conditions 
u0 = 0 
uf = 0 

#Convert to 0-based indexing
J -= 1


# In[5]:


#create arrays 
u = np.zeros(len(x))+u0
eta = np.zeros(len(x))+eta0

u_new = np.zeros(len(x))

eta_new = np.zeros(len(x))


# In[6]:


#FTCS Algorithm 
plt.figure(dpi = 150)

#start time 
t = ti
while t < tf:
    Fu,Fn = flux_cons(u,eta,H) 
    const = dt/(2*dx)
    
    #apply forward difference method for first point
    u_new[0] = u0 #Provided boundary condition
    eta_new[0] = eta[0] - 2*const*(Fn[1] - Fn[0]) 
    
    #apply central difference method for middle points
    for j in range(1,J): 
        u_new[j] = u[j] - const*(Fu[j+1] - Fu[j-1]) 
        eta_new[j] = eta[j] - const*(Fn[j+1] - Fn[j-1])
    
    #apply backwards difference method for end point
    u_new[J] = 0 #Provided boundary condition
    eta_new[J] = eta[J] - 2*const*(Fn[J] - Fn[J-1]) 
    
    #update arrays
    u = np.copy(u_new)
    eta = np.copy(eta_new)
    
    #Make plots at the given times
    if abs(t-t1) < eps:
        plt.plot(x,eta, label = f't = {t1} s')
    if abs(t-t2) < eps:
        plt.plot(x, eta, label = f't = {t2} s')
    if abs(t-t3) < eps:
        plt.plot(x, eta,label = f't = {t3} s')
        #print(eta)
     #step time
    t += dt 
    
plt.title('1D Shallow Water System with the FTCS Method')       
plt.xlabel('x')
plt.ylabel(r'$\eta$')
plt.xlim(xi,xf)
plt.legend()
plt.savefig('L08Q2.png', dpi = 150)
plt.show()


# In[ ]:




