#!/usr/bin/env python
# coding: utf-8

# In[7]:


""" Written by Madeline Nardin September 2022 """
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c,h,k, sigma
import sympy


# In[8]:


#Q3(b) Pseudocode 
#i. Define function for C1 as a function of T derived in 3(a) C1(T) = ((2*pi(kT)^4)/(c^2h^3))
#ii. Define function for integrand shown in  3(a) f(x) = (x^3/(e^x-1))
#iii. Solve integral of f(x) with Simpson method approximating bounds from x= 0.0001 to x = 709
#iv. Compute pracical error of Simpson method for N =1000 and N=10000
#iii. Solve integral of f(x) from x= 0 to x= +inf with scipy.integrate.quad() function and \
#print solution and estimated error
#iv. Multiply integral by C1(T) with desired temperature value to determine W(T)
#v. Determine value of Stefan-Boltzmann Constant by rearranging eq. 6 sigma = W*T^{-4}
#vi. Print derived sigma and true sigma (from scipy.constants)


# In[9]:


#define C1 coefficient from q3(a)
def C1(T):
    """ Determines the value of C1.

            Determines the value of C1 as a function of temperature with the following equation\
            derived in 3(a) C1(T) = ((2*pi(kT)^4)/(c^2h^3)).
        Parameters
        --------
        T : float
            T is the desired temperature.

        Returns
        -------
        C1 : float
            C1 is the constant found outside the total Energy integral derived in 3(a).
    """
    return ((2*np.pi*k**4*T**4)/(c**2*h**3))

#define integrand function 
def f(x):
    """ Defines the integrad of the total energy integral derived in 3(a) given by         f(x) = x^3/(e^x-1).
    """
    return (x**3)/((np.exp(x))-1)

def sb_const(w,T):
    """ Determine value of Stefan-Boltzmann Constant.

        Determine value of Stefan-Boltzmann Constant. by rearranging eq. 6 \
        sigma = W*T^{-4}
        Parameters
        --------
        W : float
            W is the total energy.
        
        T : float
            T is the temperature.

        Returns
        -------
        sigma : float
            sigma is the value of the Stefan-Boltzmann constant approximated by the system parameters.
    """
    sigma = w*T**(-4)
    print('numerical sigma = ', sigma, 'Wm^{-2}K^{-4}')
    return sigma


# In[56]:


def simpson(f, N, a, b):
    """ Author: Milica Ivetic
    This function uses Simpson's rule to evaluate a given function
    f from point a to b, using N slices. 
    INPUT:
    f is the function being considered in the question
    N [integer] is the number of slices 
    a, b [integer] are the starting points and end points of the integral
    
    OUTPUT:
    res [float] value of the definite integral"""    
    
    h = (b-a)/N
    
    s = f(a) + f(b)
    
    s_odd = 0
    s_even = 0
    
    for k in range(1, N, 2):
        s_odd += f(a + k*h)
        
    for k in range(2, N, 2):
        s_even += f(a +k*h)
    
    res = (1/3)*h*(s + 4*s_odd + 2*s_even)
    
    return res
#Define integratiion parameters
N1 = 10000 #number of slices 
N2 = 100000
a = 0.00001 #lower bound - note x!=0 due to zero denominator error
b = 709 #upper bound - highest value without overflow in exponent 

#checking integral with wolfram alpha's result
integral_N1 = simpson(f, N1, a, b)
print('integral:', integral_N1)

integral_N2 = simpson(f, N2, a, b)
print('integral:', integral_N2)

practical_err = (1/15)*(integral_N2 - integral_N1)
print('Practical Error = ', practical_err)


# In[57]:


from scipy.integrate import quad
#Solve integral of f(x) from x=0 to x=+inf print solution and esitmated error
I = quad(f, a=0, b=np.inf)
print('Integral = ', I[0])
print('Estimated Error = ', I[1])


# In[58]:


#Determine total energy by multiplying the C1 by the integral solution from integrate.quad()
W_500 = C1(500)*I[0]
W_800 = C1(800)*I[0]
W_1200 = C1(1200)*I[0]

#Determine total energy by multiplying the C1 by the integral solution from Simpsons Method
Ws_500 = C1(500)*integral_N2
Ws_800 = C1(800)*integral_N2
Ws_1200 = C1(1200)*integral_N2


# In[59]:


#Compute real Stefan-Boltzmann constant and approximation for T=500, T=800, and T=1200
s_500 = sb_const(W_500, 500)
s_800 = sb_const(W_800, 800)
s_1200 = sb_const(W_1200, 1200)
print('Real sigma = ', sigma)


# In[ ]:


#Note numerical sigma is accurate until the 5th decimal point with use of integrate.quad()


# In[60]:


s_500 = sb_const(Ws_500, 500)
s_800 = sb_const(Ws_800, 800)
s_1200 = sb_const(Ws_1200, 1200)
print('Real sigma = ', sigma)


# In[ ]:


#Note numerical sigma is accurate until the 5th decimal point with use of Simpsons method

