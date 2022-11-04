#This code uses the Gauss-Seidel method (with and without overrelaxation)
#to calculate the electrostatic potential at each grid point in the simple model 
# of an electronic capacitor. 
#Question 1a is modified code from Newman's (2012) laplace.py.

#Author: Milica Ivetic, November 2022

import numpy as np
import matplotlib.pyplot as plt
from pylab import imshow, gray, show
import scipy.integrate as integrate
from time import time

plt.rcParams['axes.titlesize'] = 13
plt.rcParams['axes.labelsize'] = 13
plt.rcParams['xtick.labelsize'] = 11
plt.rcParams['ytick.labelsize'] = 11


#Question 1a
print('Question 1a - see plot')

#Pseudocode for question 1a
#Make a function that does the following:
#1. Define constants, initialize arrays needed
#2. Use Gauss-Seidel method to update the gridpoints that do not have fixed voltages. 
#3. Plot (not part of the function)

def GS():
    """This function uses the Gauss-Seidel method to solve for phi in the grid."""
    
    #1. Define constants, initialize arrays needed
    M = 100 #grid squares on a side
    V_left = 1.0 #voltage of left plate
    V_right = -1.0 #voltage of right plate
    target = 1e-6 #Target accuracy

    phi = np.zeros([M+1, M+1], float)
    phi_old = np.zeros([M+1, M+1], float) #array to store old phi values
    phi[20:80, 20] = V_left #rows 20-80, 20th column has left potential 
    phi[20:80, 80] = V_right #rows 20-80, 80th column has right potential


    #Main loop
    delta = 1.0

    while delta > target:
    
        #Calculate new values of the potential
        for i in range(M+1):
            for j in range(M+1):
                if i == 0 or i == M or j == 0 or j == M: #phi is fixed
                    phi_old[i,j] = phi[i,j]
                    phi[i,j] = phi[i,j]
                
                elif i in range(20,80) and j in [20,80]: #phi is fixed
                    phi_old[i,j] = phi[i,j]
                    phi[i,j] = phi[i,j]
                
                else: #else, update phi
                    phi_old[i,j] = phi[i,j]
                    phi[i,j] = (1/4)*(phi[i+1, j] + phi[i-1, j] + phi[i, j+1] + phi[i, j-1])
                
        #Calculate maximum difference from old values
    
        delta = np.max(np.abs(phi_old - phi))
        #print(delta) , ignore
        
    return phi


#3. Contour Plot

x = np.arange(0,101)
y = np.arange(0,101)

start1 = time() #tracking how long this method takes, to compare with overrelaxation later
phi = GS()
end1 = time()

diff1 = end1 - start1

print('Gauss-Seidel method takes', diff1, 'seconds to complete for this case.')

plt.figure()
plt.contourf(x, y, phi, 30, cmap = 'twilight')
plt.colorbar(label = 'Voltage [V]')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Contour plot of potential')
plt.show()



#Question 2b
print('Question 2b - see plots')

#Pseudocode for question 2b
#1. Repeat the same steps as in part a, the only thing that changes is the 
# incorporation of the overrelaxation factor w when updating phi. 


#function for GS method with overrelaxation
def GS_overrelaxation(w):
    """ This function takes input parameter w, which changes the level of relaxation.
    It then uses Gauss-Seidel method to calculate phi for the grid."""
    
    #1. Define constants, initialize arrays needed
    M = 100 #grid squares on a side
    V_left = 1.0 #voltage of left plate
    V_right = -1.0 #voltage of right plate
    target = 1e-6 #Target accuracy

    phi = np.zeros([M+1, M+1], float)
    phi_old = np.zeros([M+1, M+1], float)
    phi[20:80, 20] = V_left #rows 20-80, 20th column has left potential 
    phi[20:80, 80] = V_right #rows 20-80, 80th column has right potential


    #Main loop
    delta = 1.0

    while delta > target:
    
        #Calculate new values of the potential
        for i in range(M+1):
            for j in range(M+1):
                if i == 0 or i == M or j == 0 or j == M:
                    phi_old[i,j] = phi[i,j]
                    phi[i,j] = phi[i,j]
                
                elif i in range(20,80) and j in [20,80]:
                    phi_old[i,j] = phi[i,j]
                    phi[i,j] = phi[i,j]
                
                else:
                    phi_old[i,j] = phi[i,j]
                    phi[i,j] = ((1+w)/4)*(phi[i+1, j] + phi[i-1, j] + phi[i, j+1] + phi[i, j-1]) - w*phi[i,j]
                
        #Calculate maximum difference from old values
    
        delta = np.max(np.abs(phi_old - phi))
        #print(delta) , ignore
        
    return phi


#Run GS_overrelaxation function for w = 0.1 and w = 0.5, plot both

#w = 0.1
start2 = time()
phi_w_a = GS_overrelaxation(0.1)
end2 = time()
diff2 = end2 - start2
print('GS-overrelaxation method takes', diff2, 'seconds to complete when w = 0.1.')

plt.figure()
plt.contourf(x, y, phi_w_a, 30, cmap = 'twilight')
plt.colorbar(label = 'Voltage [V]')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Contour plot of potential, w = 0.1')
plt.show()


#w = 0.5
start3 = time()
phi_w_b = GS_overrelaxation(0.5)
end3 = time()
diff3 = end3 - start3
print('GS-overrelaxation method takes', diff3, 'seconds to complete when w = 0.5.')

plt.figure()
plt.contourf(x, y, phi_w_b, 30, cmap = 'twilight')
plt.colorbar(label = 'Voltage [V]')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Contour plot of potential, w = 0.5')
plt.show()



#w = 0.9
start4 = time()
phi_w_b = GS_overrelaxation(0.9)
end4 = time()
diff4 = end4 - start4
print('GS-overrelaxation method takes', diff4, 'seconds to complete when w = 0.9.')

