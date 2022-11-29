#!/usr/bin/env python
# coding: utf-8

# In[1]:


''' Written by Madeline Nardin November 2022 '''
import numpy as np 
import matplotlib.pyplot as plt
from random import random


# In[2]:


#Define constants 
N = 10000 #number of sample points
T = 100 #number of computations 

a = 0 #lower limit
b = 1 #upper limit


# In[3]:


#Define integrand
def int(x):
    return x**(-1/2)/(1+np.exp(x))


# In[4]:


#write function to compute integral via mean value method n times
def compute_mean_value_mthd(T):
    sol = np.zeros(100)
    for i in range(T):
        int_vals = []
        for k in range(N):
            x = random()
            k = int(x)
            int_vals.append(k)

        sol[i] = (b-a)*np.sum(int_vals)/N
    return sol


# In[5]:


mean_method_sol = compute_mean_value_mthd(T)


# In[6]:


#Define weighting function 
def wfunc(x):
    return x**(-0.5)

#Define probability density function
def P(x):
    return 1/(2*np.sqrt(x))


# In[10]:


#Write function to compute integral via importance sampling method n times
def compute_importance_sampling_mthd(T):
    sol = np.zeros(100)
    for i in range(T):
        int_vals = []
        for k in range(N):
            x = (random())**2
            g = int(x)/wfunc(x)
            int_vals.append(g)

        sol[i] = (2/N)*np.sum(int_vals)
    return sol


# In[11]:


importance_sampling_sol = compute_importance_sampling_mthd(T)


# In[15]:


#Plot histograms of solutions from each method using 10 bins from 0.80 to 0.88
fig, ax = plt.subplots(1, 2, dpi = 150)#, sharey=True)
plt.suptitle(r'Integral of $f(x) = \frac{x^{-1/2}}{1+e^{x}}$ from $x=0$ to $x=1$')
ax[0].set_title('Mean Value Solution')
ax[0].hist(mean_method_sol, 10, range=[0.8, 0.9])
ax[0].set_xlim(0.8, 0.88)

ax[1].set_title('Importance Sampling Solution')
ax[1].hist(importance_sampling_sol, 10, range=[0.8, 0.9])
ax[1].set_xlim(0.8, 0.88)

fig.tight_layout()
fig.savefig('L10_Q3.png', dpi = 150, bbox_inches='tight')


# In[ ]:




