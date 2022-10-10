#!/usr/bin/env python
# coding: utf-8

# In[1]:


""" Written by Madeline Nardin October 2022 """
import numpy as np 
import matplotlib.pyplot as plt


# In[2]:


#Define constants 
a = 1.60218E-18 #J = 10eV
L = 5/1e+10 #m
m_e = 9.1094E-31 #Kg
q_e = 1.6002E-19 #C
hbar = 1.054571800E-34  # Js


# In[3]:


#2b
def Hmn(m,n, p):
    """ Calculates values elements of Hamiltonian matrix H according to eq.5 in assignment sheet.
    
     Parameters
        --------
        m : int
            m

        n : int 
            n
            
        p : string 
            p == 'True' includes print statements 
            p == 'False' excludes print statements

        Returns
        -------
        H : array
            H is the value of element [m,n] in the Hamiltonian matrix.
    """
    if m != n and (((n % 2) == 0 and (m % 2) == 0) or ((n % 2) != 0 and (m % 2) != 0)):
        H = 0 
        if p== 'True': print('m != n and m and n are both even or both odd')
    if m != n and (((n % 2) == 0 and (m % 2) != 0) or ((n % 2) != 0 and (m % 2) == 0)):
        H = ((8*a*m*n)/(np.pi**2*(m**2-n**2)**2))
        if p== 'True': print('m != n and one is even one is odd')
    if m == n:
        H = 0.5*a+((np.pi**2*hbar**2*m**2)/(2*m_e*L**2))
        if p== 'True': print('m = n')
    return H


# In[4]:


#2c
#define 10x10 H array
mmax = 10
nmax = 10


# In[5]:


def Hmatrix(m, n):
    """ Creates matrix of size mmax and nmax (previously defined) and stores values in python 
    style indexing.
    
     Parameters
        --------
        m : int
            m

        n : int 
            n

        Returns
        -------
        H : array
            H is the list of elements in the matrix.
    """
    #create an enpty matrix of size mmax x nmax 
    H = np. empty ([mmax, nmax], float)
    #shift index to store from 0 
    for m in range(1, mmax+1):
         for n in range(1, nmax+1):
                #calculate hamiltonian matrix elements with shifted index value
                H[m-1, n-1] = Hmn(m, n, p = 'False')
    #return matrix elements in units of eV
    return H*6.242E18 #convert Joules to eV


# In[6]:


#Compute the eigenvalues of the matrix H size 10x10
np.linalg.eigvalsh(Hmatrix(0,0))


# In[7]:


print('Ground State Energy = ', np.linalg.eigvalsh(Hmatrix(0,0))[0], 'eV')


# In[8]:


#2d
#define 100x100 H array
mmax = 100
nmax = 100


# In[9]:


#Compute the eigenvalues of the matrix H size 100x1000
np.linalg.eigvalsh(Hmatrix(0,0))[0:10]


# In[10]:


print('Ground State Energy = ', np.linalg.eigvalsh(Hmatrix(0,0))[0], 'eV')


# In[11]:


#2e
#derive eigenvalues and vectors from 100x100 matrix H
E_n, psi_n = np.linalg.eigh(Hmatrix(0,0))

def psi(level,x):
    """ Calculates the wavefunction psi(x) given the eigenvector of the matrix H from the equation 
    given in the text.
    
     Parameters
        --------
        level : int
            Energy level at which psi is computed.

        x : array 
            x is the domain in which psi si calculated over.

        Returns
        -------
        psi : array
            psi is all the points of the wavefunction over the defined domain x. 
    """
    psi = 0
    for n in range(nmax):
        psi += psi_n[level,n]*np.sin((np.pi*n*x)/(L))
    return psi


# In[12]:


def prob_density(level, x):
    return np.conjugate(psi(level, x))*psi(level, x)


# In[13]:


def normalized_psi(level):
    """ This function uses Simpson's rule to evaluate a given function
    f from point a to b, using N slices and uses the result to normalize the wavefunction psi.
    Parameters
        --------
        level : int
            Energy level at which psi is computed.

        Returns
        -------
        norm_int : array
            norm_int is all the points of the normalized wavefunction over the defined domain x. 
    """
    N = 100
    a = 0.0
    b = L
    h = (b-a)/N
    
    s = prob_density(level,a) + prob_density(level,b)
    
    s_odd = 0
    s_even = 0
    
    for k in range(1, N, 2):
        s_odd += prob_density(level, a + k*h)
        
    for k in range(2, N, 2):
        s_even += prob_density(level, a +k*h)
    
    res = (1/3)*h*(s + 4*s_odd + 2*s_even)

    print('Unnormalized Integral = ', res)
    norm_c = np.sqrt(res)
    
    norm_int = prob_density(level, x)/norm_c
    #print('Normalized Integral = ', norm_int)
    return norm_int
    


# In[14]:


#calculate and plot prob density at first three states
x = np.linspace(0,L, len(psi_n))

plt.figure(figsize = (10,8))
plt.plot(x,normalized_psi(0),'.-', label = 'Ground')
plt.plot(x,normalized_psi(1),'.-', label = 'First Excited State')
plt.plot(x,normalized_psi(2),'.-', label = 'Second Excited State')

plt.title('Normalized wavefunctions of First Three States', fontsize = 15)
plt.xlabel('x [m]', fontsize = 15)
plt.ylabel(r'$|\psi(x)|^2$', fontsize = 15)
plt.legend()

plt.xlim(0,L)
plt.savefig('L04_Q2e.png')

