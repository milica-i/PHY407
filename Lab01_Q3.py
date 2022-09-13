# Timitng Matrix multiplication
# This code ___
# Author: Milica Ivetic

import numpy as np
from time import time

# Initialize a time array

times = np.zeros(6, float)

# Using the code snippet from example 4.3 of the textbook

N_range = [2, 10, 50, 100, 200, 500]

# Iterate over range of Ns. We will time how long each iteration takes. 
For N in N_range:
    
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
    
    np.append(times, diff) 
          
print(times)      
                
# Same operation but using np.dot
    
    