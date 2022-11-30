#This code simulates Brownian motion and diffusion limited aggregation. 

#Author: Milica Ivetic, November 2022

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['axes.titlesize'] = 13
plt.rcParams['axes.labelsize'] = 13
plt.rcParams['xtick.labelsize'] = 11
plt.rcParams['ytick.labelsize'] = 11

#Question 1a, Question 10.3 pg. 457 from Newman (2012)
print('Question 1a - see plots')

#Pseudocode for Question 1a
#At each timestep, update position of the particle using nextmove function. 

#The following code is modified from Brownian-start.py provided.

"""
This program simulates Brownian motion in the presence of walls
Note that the physical behaviour would be to stick to walls,
which is the purpose of Q1a.
Author: Nico Grisouard, University of Toronto
"""

def nextmove(x, y):
    """ randomly choose a direction
    0 = up, 1 = down, 2 = left, 3 = right"""
    direction =  np.random.randint(0,4)
    if direction == 0:  # move up
        y += 1
    elif direction == 1:  # move down
        y -= 1
    elif direction == 2:  # move right
        x += 1
    elif direction == 3:  # move left
        x -= 1
    else:
        print("error: direction isn't 0-3")
    return x, y


Lp = 101  # size of domain
Nt = 5000  # number of time steps

# arrays to record the trajectory of the particle
xs = np.zeros(len(range(Nt)))
ys = np.zeros(len(range(Nt)))


centre_point = (Lp-1)//2  # middle point of domain
xp = centre_point
yp = centre_point

xs[0] = xp
ys[0] = xp
for i in range(1,Nt):
    xpp, ypp = nextmove(xp, yp)
    xp, xs[i] = xpp, xpp
    yp, ys[i] = ypp, ypp
    
#plot

plt.figure()
plt.scatter(range(Nt), xs)
plt.title('x position of particle over time')
plt.xlabel('t')                         
plt.ylabel('x')
plt.show()

plt.figure()
plt.scatter(range(Nt), ys)
plt.title('y position of particle over time')
plt.xlabel('t')
plt.ylabel('y')
plt.show()

plt.figure()
plt.scatter(xs, ys)
plt.title('x vs. y positions of particle')
plt.xlabel('x')
plt.ylabel('y')
plt.show()




#Question 1b, Question 10.13a pg.499-500 from Newman (2012)
print('Question 1b')

#Pseudocode for Question 1b
#For the specified number of particles, do the following:
#assigned new particle to start at centre of the grid
#while the centre of the grid is not an anchored point
#    if the point is at the edge, or next to an anchored point, add
#    that particle to the list of anchored points, and break the while loop
#    in order to move on to the next particle
#    otherwise, update the particle's position with nextmove function
#Plot at the end

#The following code is modified from DLA-start.py provided.

"""
This program simulates diffusion limited aggregation on an LxL grid.
Particles are initiated until the centre point is filled.
Author: Nico Grisouard, University of Toronto
Based on Paul J Kushner's DAL-eample.py
"""


Lp = 101  # size of domain
N = 1250  # number of particles


# list to represent x and y positions of anchored points
anchored_points = [] 
anchored_x = []
anchored_y = []
centre_point = (Lp-1)//2  # middle point of domain


for i in range(N):
  #assigned new particle to start at centre of the grid
    xp = centre_point
    yp = centre_point
   

    #While there are no anchored particles in the centre
    while (centre_point, centre_point) not in anchored_points:
        
        
        # if the point is at the edge, or next to an anchored point, add
        # that particle to the list of anchored points, and break the while loop
        # in order to move on to the next particle
        # otherwise, update the particle's position with nextmove function        
                        
        if xp <= 1 or xp >= Lp-1:
            anchored_x.append(xp)
            anchored_y.append(yp)
            anchored_points.append((xp, yp))
            break 
        elif yp <= 1 or yp >= Lp-1:
            anchored_x.append(xp)
            anchored_y.append(yp)
            anchored_points.append((xp, yp))    
            break
        elif (xp + 1, yp) in anchored_points:
            anchored_x.append(xp)
            anchored_y.append(yp)
            anchored_points.append((xp, yp))
            break
        elif (xp - 1, yp) in anchored_points:
            anchored_x.append(xp)
            anchored_y.append(yp)
            anchored_points.append((xp, yp))
            break            
        elif (xp,yp + 1) in anchored_points:
            anchored_x.append(xp)
            anchored_y.append(yp)
            anchored_points.append((xp, yp))
            break            
        elif (xp, yp - 1) in anchored_points:
            anchored_x.append(xp)
            anchored_y.append(yp)
            anchored_points.append((xp, yp))
            break   
        else:
            xp, yp = nextmove(xp, yp)
                

#Plot

plt.figure()
plt.scatter(anchored_x, anchored_y, color='black')
plt.title('Anchored points')
plt.xlabel('x')
plt.ylabel('y')
plt.xlim(0,101)
plt.ylim(0,101)
plt.show()