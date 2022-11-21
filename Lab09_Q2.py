#This code simulates a 2D resonant cavity being excited by imposing a z-directed 
#alternating current pattern of angular frequency w. We calculate the generated
#magnetic and electric field patterns.This is done by finding 
#the spectral solutions to Maxwell Equations. 

#Author: Milica Ivetic, November 2022

import numpy as np
import matplotlib.pyplot as plt
from dcst_for_q2 import dst, dct, idst, idct
from dcst import dst2, idst2

plt.rcParams['axes.titlesize'] = 13
plt.rcParams['axes.labelsize'] = 13
plt.rcParams['xtick.labelsize'] = 11
plt.rcParams['ytick.labelsize'] = 11


#Question 2a
print('Question 2a')
#Pseudocode for question 2a:
#1. Using dcst.py and dcst_for_q2.py, write out functions to compute 2D FT and inverses
#Only need to write functions for dcst2, dsct2 and their inverses, because the rest were imported.
#2. Test the functions that weren't given with test array f

#Function for dcst2

def dcst2(y):
    M = y.shape[0]
    N = y.shape[1]
    a = np.empty([M,N],float)
    b = np.empty([M,N],float)

    for i in range(M):
        a[i,:] = dct(y[i,:])
    for j in range(N):
        b[:,j] = dst(a[:,j])

    return b

#Function for idcst2

def idcst2(b):
    M = b.shape[0]
    N = b.shape[1]
    a = np.empty([M,N],float)
    y = np.empty([M,N],float)

    for i in range(M):
        a[i,:] = idct(b[i,:])
    for j in range(N):
        y[:,j] = idst(a[:,j])

    return y


#Function for dsct2

def dsct2(y):
    M = y.shape[0]
    N = y.shape[1]
    a = np.empty([M,N],float)
    b = np.empty([M,N],float)

    for i in range(M):
        a[i,:] = dst(y[i,:])
    for j in range(N):
        b[:,j] = dct(a[:,j])

    return b


#Function for idsct2

def idsct2(b):
    M = b.shape[0]
    N = b.shape[1]
    a = np.empty([M,N],float)
    y = np.empty([M,N],float)

    for i in range(M):
        a[i,:] = idst(b[i,:])
    for j in range(N):
        y[:,j] = idct(a[:,j])

    return y

#2. Test functions above with test array f

#Define test array

f = np.random.rand(3, 3)
f[0,:] = 0
f[:,0] = 0

print('Test array f is:', f)
print('Running the test for (i)dsct2 returns:', idsct2(dsct2(f)))
print('Running the test for (i)dcst2 returns:', idcst2(dcst2(f)))


#Question 2b
print('Question 2b - see plots')
#Pseudocode for 2b
#1. Initialize Jz, Hx, Hy, Ez = 0 everywhere at t=0 and set constants
#In a loop over time:
#2. Use eqn 7 from the lab to generate Jz(x, y, t)
#3. Take fourier transforms of Hx, Hy, Jz, and Ez
#4. Evolve Fourier coeffs using eqns. 16
#5. Reconstruct Hx, Hy, and Ez via 2D inverse FT
#6. Plot


#1. Initialize Jz, Hx, Hy, Ez = 0 everywhere at t=0 and set constants 

tau = 0.01 #timestep
N = 2000 #number of timesteps
T = N*tau #final time

#The following commented out constants don't need to be saved since they are just 1
#Lx = 1
#Ly = 1
#J0 = 1
#m = 1
#n = 1
#c = 1
P = 32
w = 3.75 #driving frequency
pi = np.pi

#discretizing t, x, and y as outlined in Algorithmic Background section
times = np.arange(0, T, tau)

ax = 1/P #size of grid cells in x-direction
ay = 1/P #size of grid cells in y-direction

p = np.arange(1, P-1)#, ax
q = np.arange(1, P-1)#, ay
x = ax*p
y = ay*q


Jz = np.zeros([len(x), len(y)], float)
Hx = np.zeros([len(x), len(y)], float)
Hy = np.zeros([len(x), len(y)], float)
Ez = np.zeros([len(x), len(y)], float)

#Fourier coefficients that get updated in the loop:
E_coeff = np.zeros([len(x),len(y)], float)
X_coeff = np.zeros([len(x),len(y)], float)
Y_coeff = np.zeros([len(x),len(y)], float)

#Keeping track of Hx, Hy and Ez traces over time
Hx_trace = []
Hy_trace = []
Ez_trace = []

#Define D, in the lab we need to use Dx and Dy but for this calculation Dx = Dy
D = pi*tau/2


#loop over time
for i in range(len(times)):
    #use eqn 7 to generate current pattern Jz(x, y, t)
    for j in range(len(x)):
        for k in range(len(y)):
            Jz[j,k] = np.sin(pi*x[j])*np.sin(pi*y[k])*np.sin(w*times[i])
    
    #Take FTs of Hx, Hy, Jz, and Ez
    Hx_FT = dsct2(Hx)
    Hy_FT = dcst2(Hy)
    Jz_FT = dst2(Jz)
    Ez_FT = dst2(Ez)
    #print(i)
    
    #Evolve the Fourier coefficients using eqns 16
    
    for m in range(len(x)):
        for n in range(len(y)):
            E_coeff[m,n] = ((1 - (m**2)*D**2 - n**2*D**2)*Ez_FT[m,n] + 2*n*D*Hx_FT[m,n] + 2*m*D*Hy_FT[m,n] + tau*Jz_FT[m,n])/(1 + (m**2)*D**2 + (n**2)*D**2)
            
            X_coeff[m,n] = Hx_FT[m,n] - n*D*(E_coeff[m,n] + Ez_FT[m,n])
            
            Y_coeff[m,n] = Hy_FT[m,n] - m*D*(E_coeff[m,n] + Ez_FT[m,n])
    
    
    
    #Reconstruct Hx, Hy, and Ez via 2D inverse FT
    
    Hx = idsct2(X_coeff)
    Hy = idcst2(Y_coeff)
    Ez = idst2(E_coeff)
    
    #Saving traces to plot over time later
    Hx_trace.append(Hx[int(len(x)/2), 1])
    Hy_trace.append(Hy[1,int(len(y)/2)])
    Ez_trace.append(Ez[int(len(x)/2), int(len(y)/2)])
    
    
    
#Plot

plt.figure()
plt.title('$E_z$(x = 0.5, y = 0.5) as a function of time')
plt.scatter(times, Ez_trace)
plt.ylabel('$E_z$')
plt.xlabel('Time [s]')
plt.show()


plt.figure()
plt.title('$H_x$(x = 0.5, y = 0.0) as a function of time')
plt.scatter(times, Hx_trace)
plt.ylabel('$H_x$')
plt.xlabel('Time [s]')
plt.show()

plt.figure()
plt.title('$H_y$(x = 0.0, y = 0.5) as a function of time')
plt.scatter(times, Hy_trace)
plt.ylabel('$H_y$')
plt.xlabel('Time [s]')
plt.show()
