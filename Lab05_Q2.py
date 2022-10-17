#!/usr/bin/env python
# coding: utf-8

# In[ ]:


""" Written by Madeline Nardin October 2022 """
""" The following code applies a “low pass filter” to the sound file 'GraviteaTime.wav'."""


# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import fft, ifft, fftfreq


# In[2]:


#Q2b Pseudocode 
#i. Read in the audio file 'GraviteaTime.wav' and separate the channels into different arrays with 
    #the code provided in the manual.
#ii. Define the sample frequency and determine the total time of the audio clip 
    #t_total = # samples/sample frequency.
#iii. Define a time array with length of the total time and same number of samples as the audio data
#iv. Plot each channel data with the time array.


# In[3]:


from scipy.io.wavfile import read, write
# read the data into two stereo channels
# sample is the sampling rate, data is the data in each channel,
#  dimensions [2, nsamples]
sample, data = read('GraviteaTime.wav')
# sample is the sampling frequency, 44100 Hz # separate into channels
channel_0 = data[:, 0]
channel_1 = data[:, 1]
N_Points = len(channel_0)


# In[4]:


#define sample rate 
sample_freq = sample
n_samples_0 = len(channel_0)
n_samples_1 = len(channel_1)
if n_samples_0 == n_samples_1: n_samples = len(channel_0)

#calculate total time of audio file in seconds
t_total = n_samples/sample_freq

#define time array 
t = np.linspace(0,t_total, n_samples)


# In[5]:


#plot channels as a function of time
fig, axs = plt.subplots(2, dpi = 150)
fig.suptitle('GraviteaTime.wav as a Function of Time')
axs[0].plot(t, channel_0)
axs[0].set_ylabel('Channel 0')
axs[0].set_xlabel('Time (s)')

axs[1].plot(t, channel_1)
axs[1].set_ylabel('Channel 1')
axs[1].set_xlabel('Time (s)')

fig.tight_layout()
fig.savefig('L05Q2b')


# In[17]:


#define a focused time interval approx 30ms long
t_focus = t[0:1500]
#plot each channel signal over this smaller interval 
fig, axs = plt.subplots(2, dpi = 150)
fig.suptitle('GraviteaTime.wav over first 30 ms')
axs[0].plot(t_focus, channel_0[0:1500])
axs[0].set_ylabel('Channel 0')
axs[0].set_xlabel('Time (s)')

axs[1].plot(t_focus, channel_1[0:1500])
axs[1].set_ylabel('Channel 1')
axs[1].set_xlabel('Time (s)')


fig.tight_layout()
fig.savefig('L05Q2c')


# In[7]:


#Q2d Pseudocode 
#i. Define frequency array with the relation ...
#ii. Fourier transform the signal to the frequency domain with np.fft().
#iii. Apply low-pass filter by setting Fourier coefficients for frequency values > 880 Hz to zero.
#iv. Take the inverse Fourier transformation of the filtered signal with np.ifft().
#v. Create quad plot (over the focused time interval in 2(c)) of lot the amplitude of the original 
    #Fourier coefficients, the amplitude of the filtered Fourier coefficients, the original time 
    #series, and the filtered time series for this interval.


# In[24]:


def lpf(signal, channel_number,figfile):
    fft_signal = fft(signal)
    N = len(fft_signal)
    n = np.arange(N)
    T = N/sample
    freq = n/T
    #freq = fftfreq(len(fft_signal), d=1/sample)

    fft_filtered = fft_signal.copy()

    for idx, f in enumerate(freq):
        if f > 800:
            fft_filtered[idx] = 0

    filtered_signal = ifft(fft_filtered)

    
    fig, axis = plt.subplots(2, 2, dpi = 150)
    fig.suptitle(f'Channel {channel_number} Plots')
    #plot FFT amplitude over focused interval 
    axis[0, 0].plot(freq, np.abs(fft_signal))
    axis[0, 0].set_ylabel('FFT apmlitude')
    axis[0, 0].set_xlabel('Frequency (Hz)')

    #plot filtered FFT amplitude over focused interval 
    axis[0, 1].plot(freq, np.abs(fft_filtered))
    axis[0, 1].set_ylabel('FFT apmlitude')
    axis[0, 1].set_xlabel('Frequency (Hz)')
    
    #plot original signal over focused interval 
    axis[1, 0].plot(t_focus, signal)
    axis[1, 0].set_ylabel('Signal Amplitude')
    axis[1, 0].set_xlabel('Time (s)')
    axis[1, 0].set_title('Original Signal')
    
    #plot filtered signal over focused interval 
    axis[1, 1].plot(t_focus, filtered_signal)
    axis[1, 1].set_ylabel('Signal Amplitude')
    axis[1, 1].set_xlabel('Time (s)')
    axis[1, 1].set_title('Filtered Signal')
    
    
    fig.tight_layout()
    fig.savefig(str(figfile))
    return filtered_signal


# In[25]:


channel_0_out = lpf(channel_0[0:1500], 0, 'L05Q2d0')


# In[26]:


channel_1_out = lpf(channel_1[0:1500], 1,'L05Q2d1')


# In[27]:


# this creates an empty array data_out with the same shape as "data" # (2 x N_Points) and the same type as "data" (int16)
data_out = np.empty((1500,2), dtype = data.dtype)
# fill data_out
data_out[:, 0] = channel_0_out
data_out[:, 1] = channel_1_out
write('GraviteaTime_lpf.wav', sample, data.astype(np.int16))


# In[ ]:




