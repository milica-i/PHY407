# Timitng Matrix multiplication
# This code follows Example 4.3 of the textbook to multiply matrices of size O(N).
# We time how long this take for various values of N. We compare these times with np.dot.
# Author: Milica Ivetic

import numpy as np
from time import time
import matplotlib.pyplot as plt

# Initialize a time array

times = []

# Using the code snippet from example 4.3 of the textbook

N_range = np.array([2, 10, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500])

# Iterate over range of Ns. We will time how long each iteration takes. 
for N in N_range:
    
#Create two constant matrices A and B using np.ones. Create array C of zeros.
    A = np.ones([N,N], float)*3 #Each entry is equal to 3
    B = np.ones([N,N], float)*4 # Each entry is equal to 4
    C = np.zeros([N,N], float)     
    
    start = time() # save start time
    
    for i in range(N): 
        for j in range(N): 
            for k in range(N): 
                C[i,j] += A[i,k]*B[k,j]
                
    end = time() # save end time
    diff = end - start # elapsed time in seconds
   
    times.append(diff)

#Plotting time as a function of N and as a function of N^3

plt.figure()
plt.plot(N_range, times)
plt.title('Matrix multiplication time as a function of N')
plt.xlabel('N')
plt.ylabel('t [sec]')
plt.show()

plt.figure()
plt.plot(N_range**3, times)
plt.title('Matrix multiplication time as a function of N^3')
plt.xlabel('N^3')
plt.ylabel('t [sec]')
plt.show()
           
                
# Same operation but using np.dot

dot_times = []

# Iterate over range of Ns. We will time how long each iteration takes. 
for N in N_range:

#Create two constant matrices A and B using np.ones. 
    A = np.ones([N,N], float)*3 #Each entry is equal to 3
    B = np.ones([N,N], float)*4 # Each entry is equal to 4  
    
    start = time() # save start time
    
    np.dot(A, B)
                
    end = time() # save end time
    diff = end - start # elapsed time in seconds
   
    dot_times.append(diff)
    
print(times, dot_times)

#Plot for comparison 

plt.figure()
plt.plot(N_range, dot_times)
plt.title('np.dot times as a function of N')
plt.xlabel('N')
plt.ylabel('t [sec]')
plt.show()

plt.figure()
plt.plot(N_range**3, dot_times)
plt.title('np.dot times as a function of N^3')
plt.xlabel('N^3')
plt.ylabel('t [sec]')
plt.show()