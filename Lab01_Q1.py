#Modelling a planetary orbit 
# This code models Mercury's orbit using the Euler-Cromer method. 
# We model consider Newtonian gravity force and gravitational force in the general relativity form. We also check to see if angular momentum is conserved. 
# Mercury's orbit is plotted (x vs. y) as well as velocities in the x and y directions over time. 
#Author: Milica Ivetic
# 

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as spc

# Question 1b (Pseudocode):

#1. Define constants: G, dt, Ms, Mp. Note we are working in AU, Ms, and Earth year units.

#2. Initialize time array for one year. Note that dt is given to us.

#3. Initialize x, y, vx, and vy arrays using given initial conditions.
# Also Initialize angular momentum arrays in both the x and y directions. Using the formula L = M*v*r where M is the mass of the planet and r is the distance 
# between the two objects which is given by (x^2 + y^2)^0.5. 

#4. For 1 year:
#       Increment the x, y, vx, and vy arrays using equations derived in 1a
#   The equations are: vx[i+1] = vx[i] - dt*(G*Ms*x[i])/(r**3)
#                      vy[i+1] = vy[i] - dt*(G*Ms*y[i])/(r**3)
#                      x[i+1] = x[i] + dt*vx[i+1]
#                      y[i+1] = y[i] + dt*vy[i+1]
# This is the Euler-Cromer method^. 

# Add the Lx and Ly components in quadrature to get an L array.

#5. Plot x vs. y, t vs. vx and t vs. vy. Check if L is conserved by plotting t vs. L



# Question 1c (Real code):

#1. Define constants:

G = 39.5 # AU^3/Ms/yr^2
Ms = 1 # Since we are working in solar masses
Mp = 1.651e-7 # Solar masses
dt = 0.0001 #yr        
        
#2. Initialize time array for one year. 

t = np.arange(0, 1, dt)

#3. Initialize x, y arrays and vx and vy arrays using initial conditions
#Make an array of zeros for x, y, vx, vy, Lx, and Ly. len(t) is the number of timesteps.

x = np.zeros(len(t))
y = np.zeros(len(t))
vx = np.zeros(len(t))
vy = np.zeros(len(t))
Lx = np.zeros(len(t))
Ly = np.zeros(len(t))

# Use initial conditions given to us in the question:

x[0] = 0.47 #AU
y[0] = 0.0 #AU, this is redundant since all values of y are zero initially
vx[0] = 0.0 # AU/yr, redundant again
vy[0] = 8.17 # AU/yr
Lx[0] = Mp*vx[0]*(np.sqrt(x[0]**2+y[0]**2)) #Using formula L = mvr
Ly[0] = Mp*vy[0]*(np.sqrt(x[0]**2+y[0]**2)) #Units are [Ms AU^2/yr]

#4. Iterate using Euler-Cromer method:

for i in range(len(t)-1):
        vx[i+1] = vx[i] - dt*(G*Ms*x[i])/(np.sqrt(x[i]**2+y[i]**2)**3)
        vy[i+1] = vy[i] - dt*(G*Ms*y[i])/(np.sqrt(x[i]**2+y[i]**2)**3)
        x[i+1] = x[i] + dt*vx[i+1]
        y[i+1] = y[i] + dt*vy[i+1]
        Lx[i+1] = Mp * (np.sqrt(x[i]**2+y[i]**2))*(vx[i+1])
        Ly[i+1] = Mp * (np.sqrt(x[i]**2+y[i]**2))*(vy[i+1])
        

# Add components of angular momentum in quadrature 

L = np.sqrt(Lx**2 + Ly**2) #[Ms AU^2/yr]

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

# Checking if angular momentum is conserved

plt.figure()
plt.plot(t, L)
plt.title('Mercury angular momentum as a function of time')
plt.xlabel('t [yr]')
plt.ylabel('L [Ms AU^2/yr]')
plt.ylim(0,0.000001)
plt.show()


plt.figure()
plt.plot(t, L)
plt.title('Mercury angular momentum as a function of time, zoomed in')
plt.xlabel('t [yr]')
plt.ylabel('L [Ms AU^2/yr]')
plt.show()



#Question 1d

#Define alpha 

a = 0.01 # AU^2

#new gravitational force law equation is :
# Fg = -GMsMp * (1 + a/r^2) * (xxhat + yyhat)


#initializeing time array
dtd = 0.0001
td = np.arange(0, 10, dtd)

#initialize x, y, vx, and vy arrays
xd = np.zeros(len(td))
yd = np.zeros(len(td))
vxd = np.zeros(len(td))
vyd = np.zeros(len(td))

#Using same initial conditions
xd[0] = 0.47 #AU
yd[0] = 0.0 #AU, this is redundant since all values of y are zero initially
vxd[0] = 0.0 # AU/yr, redundant again
vyd[0] = 8.17 # AU/yr

#Iterate, The new equations for x, y, vx, and vy are 

for i in range(len(td)-1):
        vxd[i+1] = vxd[i] - dtd*((G*Ms*xd[i])/((np.sqrt(xd[i]**2+yd[i]**2)**3))*(1 + a/(np.sqrt(xd[i]**2 + yd[i]**2)**2)))
        vyd[i+1] = vyd[i] - dtd*((G*Ms*yd[i])/((np.sqrt(xd[i]**2+yd[i]**2)**3))*(1 + a/(np.sqrt(xd[i]**2 + yd[i]**2)**2)))
        xd[i+1] = xd[i] + dtd*vxd[i+1]
        yd[i+1] = yd[i] + dtd*vyd[i+1]

#New plot of x vs. y

plt.figure()
plt.plot(xd, yd)
plt.title('x vs. y of Mercurys orbit over 10 Earth years with precession')
plt.xlabel('x [AU]')
plt.ylabel('y [AU]')
plt.show()
