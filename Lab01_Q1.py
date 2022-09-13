#Modelling a planetary orbit 
# This code ___
#Author: Milica Ivetic
# 

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as spc

# Question 1b (Pseudocode):

#1. Define constants: G, dt, r note we are working in AU, Ms, and Earth year units

#2. Initialize time array for one year. Note that dt is given to us.

#3. Initialize x, y arrays and vx and vy arrays using initial conditions

#4. For 1001 iterations:
#       Increment the x, y, vx, and vy arrays using equations derived in 1a
#   The equations are: vx[i+1] = vx[i] - dt*(G*Ms*x[i])/(r**3)
#                      vy[i+1] = vy[i] - dt*(G*Ms*y[i])/(r**3)
#                      x[i+1] = x[i] + dt*vx[i+1]
#                      y[i+1] = y[i] + dt*vy[i+1]
# This is the Euler-Cromer method^

#5. Plot x vs. y, t vs. vx and t vs. vy. 

# Question 1c (Real code):

#1. Define constants:

G = 39.5 # AU^3/Ms/yr^2
Ms = 1 # Since we are working in solar masses
dt = 0.0001 #yr
r = 0.4 #AU, this is the distance between the Sun and Mercury in AU. Acquired 
        #from online source https://www.nasa.gov/audience/foreducators/5-8/
        #features/F_Solar_System_Scale.html
        
        
#2. Initialize time array for one year. 

t = np.arange(0, 1, dt)

#3. Initialize x, y arrays and vx and vy arrays using initial conditions
#Make an array of zeros for x, y, vx, and vy .len(t) is the number of timesteps.

x = np.zeros(len(t))
y = np.zeros(len(t))
vx = np.zeros(len(t))
vy = np.zeros(len(t))

# Use initial conditions given to us in the question:

x[0] = 0.47 #AU
y[0] = 0.0 #AU, this is redundant since all values of y are zero initially
vx[0] = 0.0 # AU/yr, redundant again
vy[0] = 8.17 # AU/yr

#4. Iterate using Euler-Cromer method:

for i in range(len(t)-1):
        vx[i+1] = vx[i] - dt*(G*Ms*x[i])/(r**3)
        vy[i+1] = vy[i] - dt*(G*Ms*y[i])/(r**3)
        x[i+1] = x[i] + dt*vx[i+1]
        y[i+1] = y[i] + dt*vy[i+1]
        
        
#5. Plots

#x vs. y

plt.figure()
plt.plot(x, y)
plt.title('x vs. y of Mercurys orbit over 1 Earth year')
plt.xlabel('x [AU]')
plt.ylabel('y [AU]')
plt.show()

#t vs. vx

plt.figure()
plt.plot(t, vx)
plt.title('x-component of Mercury velocity as a function of time')
plt.xlabel('t [yr]')
plt.ylabel('vx [AU/yr]')
plt.show()

#t vs. vy

plt.figure()
plt.plot(t, vy)
plt.title('y-component of Mercury velocity as a function of time')
plt.xlabel('t [yr]')
plt.ylabel('vy [AU/yr]')
plt.show()


#Question 1d

#Define alpha 

a = 1 # AU^2

#new gravitational force law equation is :
# Fg = -GMsMp * (1 + a/r^2) * (xxhat + yyhat)

# The new equations for x, y, vx, and vy are (using same initical conditions):

dtd = 0.01
td = np.arange(0, 1000, dtd)

xd = np.zeros(len(td))
yd = np.zeros(len(td))
vxd = np.zeros(len(td))
vyd = np.zeros(len(td))

xd[0] = 0.47 #AU
yd[0] = 0.0 #AU, this is redundant since all values of y are zero initially
vxd[0] = 0.0 # AU/yr, redundant again
vyd[0] = 8.17 # AU/yr

#Iterate

for i in range(len(td)-1):
        vxd[i+1] = vxd[i] - dtd*((G*Ms*xd[i])/(r**3))*(1 + a/(r**2))
        vyd[i+1] = vyd[i] - dtd*((G*Ms*yd[i])/(r**3))*(1 + a/(r**2))
        xd[i+1] = xd[i] + dtd*vxd[i+1]
        yd[i+1] = yd[i] + dtd*vyd[i+1]

#New plot of x vs. y

plt.figure()
plt.plot(xd, yd)
plt.title('x vs. y of Mercurys orbit over 1 Earth year with precession')
plt.xlabel('x [AU]')
plt.ylabel('y [AU]')
plt.show()


        
        

        







