#This code calculates the volume of a 10-dimensional unit hypershpere using Monte Carlo integration. 

#Author: Milica Ivetic, November 2022


import numpy as np


#Question 2, exercise 10.7 on pages 471-472 from Newman (2012)

#Pseudocode for Question 2
#1. Set constants, initialize values to be updated in loop
#2. For 1 million steps, generate 10 random numbers between 0 and 1
#3. Generate 10D equation for a sphere, using generated values, take square root
# to get distance fr
#4. If square root of the above is less than 1, increment count
#5. Calculate volume using V*count/N

#1. Set constants, initialize values to be updated in loop
N = int(1e6) #number of points
d = 10 #Number of dimensions
V = 2**d #volume V in Eqn 10.33 from Newman
count = 0 #initializing value in the sum

#2. For 1 million steps, generate 10 random numbers between 0 and 1
for i in range(N):
    x = np.random.uniform(0,1,d)
    
    #initialize squared equation of 10D sphere
    R_squared = 0
    
    #3.Generate 10D equation for a sphere, using generated values, take square root
    
    for j in range(d):
        R_squared = R_squared + x[j]**2
        
    R = np.sqrt(R_squared)
    
    #4. If square root of the above is less than or equal to 1, increment count
    if R <= 1:
        count += 1


#5. Calculate volume using V*count/N

I = V*count/N #value of the integral

print('The value of the integral is:', I)