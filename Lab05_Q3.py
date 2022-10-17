#!/usr/bin/env python
# coding: utf-8

# In[8]:


""" Written by Madeline Nardin October 2022 """
""" The following code extracts component of SLP corresponding to Fourier wavenumber m by taking
    all rows of SLP fft at column m and plots in the time-longitude domain. """
import matplotlib.pyplot as plt
import numpy as np
from numpy.fft import rfft2 , irfft2


# In[2]:


SLP = np.loadtxt('SLP.txt')
Longitude = np.loadtxt('lon.txt')
Times = np.loadtxt('times.txt')


# In[3]:


plt.figure(dpi=200)
plt.contourf(Longitude, Times, SLP)
plt.xlabel('longitude(degrees)')
plt.ylabel('days since Jan. 1 2015')
plt.title('SLP anomaly (hPa)')
plt.colorbar()


# In[4]:


#Q3 Pseudocode 
#i. Take fft of SLP with np.rfft.
#ii. Define new array for extracted data.
#iii. Extract component of SLP corresponding to Fourier wavenumber m by taking
    # all rows of SLP fft at column m. 
#iv. Take inverse FFT of extracted data.
#v. Create contour plot of extracted data (colour bar) in the time-longitude domain. 


# In[5]:


def filter(m):
    #Take fft of SLP
    SLP_fft = rfft2(SLP)
    #Define array for extracted wavelength SLP
    extracted_fft= np.empty(SLP_fft.shape)
    #extract component of SLP corresponding to Fourier wavenumber m
    extracted_fft[:,m] = SLP_fft[:,m]
    #Take inverse FFT of extracted data 
    filtered_SLP = irfft2(extracted_fft)
    
    #Plot extracted data in the time-longitude domain 
    plt.figure(dpi=200)
    plt.contourf(Longitude, Times, filtered_SLP)
    plt.xlabel('Longitude [degrees]')
    plt.ylabel('Days Since Jan. 1 2015')
    plt.title(f'SLP (hPa) m = {m} ')
    plt.colorbar()
    plt.savefig(f'Q3m{m}.png')
    return filtered_SLP


# In[6]:


filter(3)


# In[7]:


filter(5)


# In[ ]:




