#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve


# In[2]:


#Define Constants
m = 9.109e-31                  
L = 1e-8                       
x0 = L/5                       
sigma = L/25
k = 500/L
hbar = 6.626*10**(-34)    

Nt = 3000                        
Nx = 1024 

dx = L/Nx                      
dt = 10**(-18) 
T = Nt*dt

x = np.arange(-L/2,L/2, dx)
t = np.arange(0,T, dt)


# In[3]:


#Define psi(x, t=0)
def psi0(x):
    a = -(x-x0)**2/(4*sigma**2)
    b = 1j*k*x
    return np.exp(a+b)


# In[4]:


#Define Potential within spatial interval
V = 0


# In[5]:


#Define psi(x)
psi = np.array(list(map(psi0, x)), complex)     # Wave equation at t = 0
#Normalize psi(x)
norm_const = np.sqrt(np.matmul(np.conj(psi), psi))   # Normalization constant
psi_normalized = psi / norm_const 


# In[6]:


#Define the Hamiltonain Matrix 
A = -hbar**2/(2*m*dx**2)
B = V + (-2*A)

#HD = np.zeros([Nx,Nx], complex)
vec_diag = np.ones(Nx)*B
D = np.diag(vec_diag, k=0)
Sup = A*np.eye(Nx, k=1)
Sub = np.conj(A*np.eye(Nx, k=-1))
HD = D + Sub + Sup; print('\n', HD)


# In[7]:


#Define Identity matrix
I= np.diag(np.ones(Nx), k=0)


# In[8]:


#Determine matrix L=(I_{P−1}+idt/2hbar HD)
Lmatrix = np.zeros([Nx,Nx], complex)
Rmatrix = np.zeros([Nx,Nx], complex)
for i in range(len(Lmatrix)):
    for j in range(len(Lmatrix)):
            Lmatrix[i][j] = I[i][j]+1j*dt*HD[i][j]/(2*hbar)
            Rmatrix[i][j] = I[i][j]-1j*dt*HD[i][j]/(2*hbar)
            
print('\n', Lmatrix)
print('\n', Rmatrix)


# In[9]:


#Loop C-N Method throguh time steps
sol = []
sol.append(psi_normalized)

# Calculate the wave equation at each time step
for i in range(Nt):
    v = np.matmul(Rmatrix, sol[-1])
    psi_sol = solve(Lmatrix, v)
    norm_const = np.sqrt(np.matmul(np.conj(psi_sol), psi_sol))
    psi_norm_sol = psi_sol / norm_const
    sol.append(psi_norm_sol)
    sol.append(psi_sol)


# In[15]:


fig, ax = plt.subplots(5, figsize = (10,8))
fig.suptitle('Time-Dependent Schrodinger Equation with the C-N Method', fontsize = 15)
fig.tight_layout()

ax[0].plot(x,sol[0].real)
ax[0].set_title(f't = 0.0 s')

ax[1].plot(x,sol[749].real)
ax[1].set_title(f't = T/4 s')

b = Nt/2
ax[2].plot(x,sol[1499].real)
ax[2].set_title(f't = T/2 s')

c = 3*Nt/4
ax[3].plot(x,sol[2249].real)
ax[3].set_title(f't = 3T/4 s')

ax[4].plot(x,sol[-1].real)
ax[4].set_title(f't = T s')


ax[4].set_xlabel(r'x [m]', fontsize = 15)
ax[2].set_ylabel(r'$\psi(x,t)$', fontsize = 15)
fig.savefig('L09_Q1.png', bbox_inches='tight',dpi = 150)


# In[11]:


pos_vals = []
energy_vals = []
norm_vals = []
t_vals = []

for i in range(Nt):
    sum_pos = np.matmul(x, np.conjugate(sol[i])*(sol[i]))  
    avg_pos = sum_pos / len(x)
    
    Hpsi = np.matmul(HD, sol[i])
    E = np.matmul(np.conjugate(sol[i]), Hpsi)
    
    norm = np.matmul(np.conjugate(sol[i]), sol[i])
    
    pos_vals.append(avg_pos)
    energy_vals.append(E)
    norm_vals.append(norm)
    
    tpoint = dt*i
    t_vals.append(tpoint)


# In[16]:


plt.figure(dpi = 150)
plt.plot(t_vals, pos_vals)
plt.title('Expectation Value of Position as a Function of Time')
plt.xlabel('t [s]')
plt.ylabel(r'<x>(t)')
plt.tight_layout()
plt.savefig('L09_Q1b.png', bbox_inches='tight', dpi = 150)


# In[17]:


plt.figure(dpi = 150)
plt.plot(t_vals, energy_vals, label = r'E(t)')
plt.plot(t_vals, norm_vals,  label = 'Normalization')
plt.xlabel('t [s]')
plt.ylabel(r'E(t)')
plt.title('Energy and Normalization Parameter of Wavefunction as a Function of Time')
plt.legend()
plt.xlim(0,T)
plt.tight_layout()
plt.savefig('L09_Q1c.png', bbox_inches='tight', dpi = 150)


# In[ ]:




