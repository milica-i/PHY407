#!/usr/bin/env python
# coding: utf-8

# In[1]:


from math import sqrt,exp
import numpy as np
from random import random,randrange, seed
import time
import matplotlib.pyplot as plt
#from visual import sphere,curve,display


# In[2]:


#tau = 1e2
#tau = 1e4
tau = 1e5


# In[3]:


Dvals = []

start = time.time()
for n in range(5):
    N = 25
    R = 0.02
    Tmax = 10.0
    Tmin = 1e-3

    # Function to calculate the magnitude of a vector
    def mag(x):
        return sqrt(x[0]**2+x[1]**2)

    # Function to calculate the total length of the tour
    def distance():
        s = 0.0
        for i in range(N):
            s += mag(r[i+1]-r[i])
        return s

    seed(10)
    # Choose N city locations and calculate the initial distance
    r = np.empty([N+1,2],float)
    for i in range(N):
        r[i,0] = random()
        r[i,1] = random()
    r[N] = r[0]
    D = distance()

    # Set up the graphics
    #display(center=[0.5,0.5])
    #for i in range(N):
    #    sphere(pos=r[i],radius=R)
    #l = curve(pos=r,radius=R/2)

    # Main loop
    t = 0
    T = Tmax

    seed(n)
    while T>Tmin:

        # Cooling
        t += 1
        T = Tmax*exp(-t/tau)

        # Update the visualization every 100 moves
        #if t%100==0:
            #l.pos = r

        # Choose two cities to swap and make sure they are distinct
        i,j = randrange(1,N),randrange(1,N)
        while i==j:
            i,j = randrange(1,N),randrange(1,N)

        # Swap them and calculate the change in distance
        oldD = D
        r[i,0],r[j,0] = r[j,0],r[i,0]
        r[i,1],r[j,1] = r[j,1],r[i,1]
        D = distance()
        deltaD = D - oldD


        # If the move is rejected, swap them back again
        if random()>exp(-deltaD/T):
            r[i,0],r[j,0] = r[j,0],r[i,0]
            r[i,1],r[j,1] = r[j,1],r[i,1]
            D = oldD
            
    Dvals.append(deltaD)
    
    xvals = []
    yvals = []

    for item in r:
        x,y = item
        xvals.append(x)
        yvals.append(y)
    
    plt.figure(dpi = 100)
    plt.plot(xvals,yvals, '.-')
    plt.title(f'Salesman Solution: tau ={tau}, seed = {n}, D = {D}')
    plt.grid()
    plt.savefig(f'L11q1_tau{tau}_seed{n}.png', dpi = 100)
    plt.show()

end = time.time()
print(end - start)


# #### 
