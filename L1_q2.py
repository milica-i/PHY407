#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as spc
""" Written by Madeline Nardin September 2022 ""

# Question 2a (Pseudocode):

#i. Define constants in base units of AU, Ms, and Earth year units

#ii. Define time array year. Note that dt is given to us.

#iii. Initialize position (x and y) arrays and velocity (vx and vy) arrays for Jupiter

#iv. Define Jupiter initial conditions

#iv. Iterate eqns derived in 1a with the Euler-Cromer method

#v. Initialize position (x and y) arrays and velocity (vx and vy) arrays for Earth

#vi. Define Earth initial conditions

#vii. Iterate eqns derived in 1a with the Euler-Cromer method

#viii. Plot Jupiter position and Earth Position.

#Define constants:
G = 39.5 # AU^3/Ms/yr^2
Ms = 1 # Since we are working in solar masses
MJ = 10E-3 #Ms
dt = 0.0001 #yr


# In[2]:


def Fg(x,y, Ms=Ms):
    """ Determines the Newtonian gravitational force keeping an object in orbit.

            Determines the Newtonian gravitational force keeping a planet in orbit with the following \
            relation F_g = (-GMsMp/r^3)(x+y).
        Parameters
        --------
        x : array
            x is the x component of the object's position in orbit.

        y : array 
            y is the x component of the object's position in orbit.

        Ms : int
            Ms is the mass of the object acting on the orbiting object. Equal to one solar mass \
            unless specified otherwise. 

        Returns
        -------
        Fg_x : float
            x-component of the  Newtonian gravitational force keeping an object in orbit.
        
        Fg_y : float
            y-component of the  Newtonian gravitational force keeping an object in orbit.
    """
    G = 39.5 # AU^3/Ms/yr^2
    r = np.sqrt(x**2 + y**2)
    Fg_x = ((-G*Ms)/(r**3))*x
    Fg_y = ((-G*Ms)/(r**3))*y
    return Fg_x, Fg_y

def orbit_sim_Jupiter(t):
    """ Determines orbit of Jupiter in the Sun Jupiter system over a given time interval.

            Determines the orbit of Jupiter in the Sun Jupiter system over a given time interval by \
            integrating equations from q1a to calculate the position and velocity of Jupiter as a \
            function of time under the Newtonian gravity force with the Euler method.
        Parameters
        --------
        t : array
            t is the total time that the object will orbit for. ex. t = np.arange(0, 10, dt) is an \
            orbit time of 10 Earth years.

        Returns
        -------
        xJ : array
            xJ is the x-component of the object's position in orbit.
        
        yJ : array
            yJ is the y-component of the object's position in orbit.
    """
    #Make an array of zeros for Jupiter's position and velocity components with an index number \
        #equal to the number of timesteps.
    xJ = np.zeros(len(t))
    yJ = np.zeros(len(t))
    vxJ = np.zeros(len(t))
    vyJ = np.zeros(len(t))
    axJ = np.zeros(len(t))
    ayJ = np.zeros(len(t))

    #Define Jupiter inital conditions provided in the question.
    xJ[0] = 5.2 # AU
    yJ[0] = 0 # AU
    vxJ[0] = 0 # AU/ Year
    vyJ[0] = 2.63 # AU/ Year

    #Iterate using Euler-Cromer method:
    for i in range(len(t)-1):
        #Determine acceleration (due to gravity) components
        axJ[i+1], ayJ[i+1] = Fg(xJ[i], yJ[i])

        #Determine velocity components
        vxJ[i+1] = vxJ[i] + dt*axJ[i]
        vyJ[i+1] = vyJ[i] + dt*ayJ[i]
        
        #Determinie position components
        xJ[i+1] = xJ[i] + dt*vxJ[i]
        yJ[i+1] = yJ[i] + dt*vyJ[i]
    return xJ, yJ


# In[3]:


#q2(a)
#Define time array for 10 Earth years 
t = np.arange(0, 10, dt) 

#determine and plot Jupiter's orbit wrt. the Sun 
plt.figure(figsize= (10,8))
plt.plot(0,0,'*y', label = 'Sun')
plt.plot(orbit_sim_Jupiter(t)[0], orbit_sim_Jupiter(t)[1], '--b', label = 'Jupiter')
plt.plot(orbit_sim_Jupiter(t)[0][0], orbit_sim_Jupiter(t)[1][0], 'ob')
plt.title('Jupiter Orbit Over 10 Earth Years')
plt.xlabel('x [AU]')
plt.ylabel('y [AU]')
plt.axis('square')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()


# In[ ]:


def three_body_sim(t, objects, MJ = 10E-3):
    """ Determines orbits of a three body system over a given time interval.

            Determines the orbit of Jupiter in the Sun Jupiter system over a given time interval by \
            integrating equations from q1a to calculate the position and velocity of Jupiter as a \
            function of time under the Newtonian gravity force with the Euler-Cromer method.
        Parameters
        --------
        t : array
            t is the total time that the object will orbit for. ex. t = np.arange(0, 10, dt) is an \
            orbit time of 10 Earth years.
            
        objects : list
            objects specifies what two objects make up the three body system (the sun is assumed to always \
            be at the origin). 
        
        MJ : float 
            MJ is the mass of Jupiter equal to 10E-3 solar mass \
            unless specified otherwise. 

        Returns
        -------
        if objects == 'Jupiter, Earth'
        xJ : array
            xJ is the x-component of Jupiter's position in orbit.
        
        yJ : array
            yJ is the y-component of the Jupiters's position in orbit.
        
        xE : array
            xE is the x-component of the Earth's position in orbit.
        
        yE : array
            yE is the y-component of the Earth's position in orbit.
            
        if objects == 'Jupiter, Asteroid'
        xJ : array
            xJ is the x-component of Jupiter's position in orbit.
        
        yJ : array
            yJ is the y-component of the Jupiters's position in orbit.
        
        xA : array
            xA is the x-component of the Asteroid's position in orbit.
        
        yA : array
            yA is the y-component of the Asteroid's position in orbit.
    """
    if objects == 'Jupiter, Earth':
        #Make an array of zeros for Jupiter's position and velocity components with an index number \
        #equal to the number of timesteps.
        xJ = np.zeros(len(t))
        yJ = np.zeros(len(t))
        vxJ = np.zeros(len(t))
        vyJ = np.zeros(len(t))
        axJ = np.zeros(len(t))
        ayJ = np.zeros(len(t))

        #Define Jupiter inital conditions provided in the question.
        xJ[0] = 5.2 # AU
        yJ[0] = 0 # AU
        vxJ[0] = 0 # AU/ Year
        vyJ[0] = 2.63 # AU/ Year

        #Iterate using Euler-Cromer method:
        for i in range(len(t)-1):

            #Determine acceleration (due to gravity) components
            axJ[i+1], ayJ[i+1] = Fg(xJ[i], yJ[i])

            #Determine velocity components
            vxJ[i+1] = vxJ[i] + dt*axJ[i]
            vyJ[i+1] = vyJ[i] + dt*ayJ[i]
            
            #Determine position components 
            xJ[i+1] = xJ[i] + dt*vxJ[i]
            yJ[i+1] = yJ[i] + dt*vyJ[i]
        
        #Make an array of zeros for Earth's position and velocity components with an index number \
        #equal to the number of timesteps.
        xE = np.zeros(len(t))
        yE = np.zeros(len(t))
        vxE = np.zeros(len(t))
        vyE = np.zeros(len(t))
        axE = np.zeros(len(t))
        ayE = np.zeros(len(t))

        axJE = np.zeros(len(t))
        ayJE = np.zeros(len(t))

        #Define Earth inital conditions provided in the question.
        xE[0] = 1.0 # AU
        yE[0] = 0 # AU
        vxE[0] = 0 # AU/ Year
        vyE[0] = 6.18 # AU/ Year

        #Iterate using Euler-Cromer method:
        for i in range(len(t)-1):

            #Determine acceleration (due to gravity) components 
            axE[i+1], ayE[i+1] = Fg(xE[i], yE[i])
            axJE[i+1], ayJE[i+1] = Fg(xE[i], yE[i],Ms=MJ)
            
            #Ensure Earth is between Sun and Jupiter 
            if xJ[i]>xE[i]:
                axE[i] = np.abs(axE[i])
            else:axE[i] = (-1)*np.abs(axE[i])

            if yJ[i]>yE[i]:
                ayE[i] = np.abs(ayE[i])
            else:ayE[i] = (-1)*np.abs(ayE[i])
            
            #Determine sum of acelleration due to Sun and acelleration due to Jupiter
            ax_tot = axE[i+1] + axJE[i+1]
            ay_tot =  ayE[i+1] +ayJE[i+1]

            #Determine velocity components 
            vxE[i+1] = vxE[i] + dt*ax_tot
            vyE[i+1] = vyE[i] + dt*ay_tot
            
            #Determine position components 
            xE[i+1] = xE[i] + dt*vxE[i]
            yE[i+1] = yE[i] + dt*vyE[i]
    
        return xJ,yJ, xE, yE
    
    if objects == 'Jupiter, Asteroid':
        #Make an array of zeros for Jupiter's position and velocity components with an index number \
        #equal to the number of timesteps.
        xJ = np.zeros(len(t))
        yJ = np.zeros(len(t))
        vxJ = np.zeros(len(t))
        vyJ = np.zeros(len(t))
        axJ = np.zeros(len(t))
        ayJ = np.zeros(len(t))

        #Define Jupiter inital conditions provided in the question.
        xJ[0] = 5.2 # AU
        yJ[0] = 0 # AU
        vxJ[0] = 0 # AU/ Year
        vyJ[0] = 2.63 # AU/ Year

        #Iterate using Euler-Cromer method:
        for i in range(len(t)-1):

            #Determine acceleration (due to gravity) comonents 
            axJ[i+1], ayJ[i+1] = Fg(xJ[i], yJ[i])

            #Determine velocity components
            vxJ[i+1] = vxJ[i] + dt*axJ[i]
            vyJ[i+1] = vyJ[i] + dt*ayJ[i]
            
            #Determine position components
            xJ[i+1] = xJ[i] + dt*vxJ[i]
            yJ[i+1] = yJ[i] + dt*vyJ[i]
        
        
        #Make an array of zeros for the Asteroid's position and velocity components with an index number \
        #equal to the number of timesteps.
        xA = np.zeros(len(t))
        yA = np.zeros(len(t))
        vxA = np.zeros(len(t))
        vyA = np.zeros(len(t))
        axA = np.zeros(len(t))
        ayA = np.zeros(len(t))

        axJA = np.zeros(len(t))
        ayJA = np.zeros(len(t))

        #Define Asteroid inital conditions provided in the question.
        xA[0] = 3.3 # AU
        yA[0] = 0 # AU
        vxA[0] = 0 # AU/ Year
        vyA[0] = 3.46 # AU/ Year

        #Iterate using Euler-Cromer method:
        for i in range(len(t)-1):

            #Determine acceleration (due to gravity) components
            axA[i+1], ayA[i+1] = Fg(xA[i], yA[i])
            axJA[i+1], ayJA[i+1] = Fg(xA[i], yA[i],Ms=MJ)
            
            #Ensure Earth is between Sun and Jupiter
            if xJ[i]>xA[i]:
                axA[i] = np.abs(axA[i])
            else:axA[i] = (-1)*np.abs(axA[i])

            if yJ[i]>yA[i]:
                ayA[i] = np.abs(ayA[i])
            else:ayA[i] = (-1)*np.abs(ayA[i])

            #Determine sum of acelleration due to Sun and acelleration due to Jupiter
            ax_tot = axA[i+1] + axJA[i+1]
            ay_tot =  ayA[i+1] + ayJA[i+1]

            #Determine velocity components
            vxA[i+1] = vxA[i] + dt*ax_tot
            vyA[i+1] = vyA[i] + dt*ay_tot
            
            #Determine position components
            xA[i+1] = xA[i] + dt*vxA[i]
            yA[i+1] = yA[i] + dt*vyA[i]
        return xJ,yJ, xA,yA


# In[ ]:


#Determine and plot orbits of Jupiter Earth system wrt. the sun over 10 Earth years
plt.figure(figsize= (10,8))
plt.plot(0,0,'*y', label = 'Sun')
plt.plot(three_body_sim(t, 'Jupiter, Earth', MJ = 10E-3)[0],three_body_sim(t, 'Jupiter, Earth', MJ = 10E-3)[1]         ,'--k', label = 'Jupiter')
plt.plot(three_body_sim(t, 'Jupiter, Earth', MJ = 10E-3)[0][0],three_body_sim(t, 'Jupiter, Earth', MJ = 10E-3)[1][0]         ,'ok')
plt.plot(three_body_sim(t, 'Jupiter, Earth', MJ = 10E-3)[2],three_body_sim(t, 'Jupiter, Earth', MJ = 10E-3)[3]         ,'--b', label = 'Earth')
plt.plot(three_body_sim(t, 'Jupiter, Earth', MJ = 10E-3)[2][0],three_body_sim(t, 'Jupiter, Earth', MJ = 10E-3)[3][0]         ,'ob')
plt.title('Orbits Over 10 Earth Years')
plt.xlabel('x [AU]')
plt.ylabel('y [AU]')
plt.axis('square')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig('Earth_Jupiter_10yrs.png')
plt.show()


# In[ ]:


#q2(b) Determine and plot orbits of Jupiter (increased mass MJ = 1000*10E-3 Ms) Earth system wrt. the sun over 3 & 3.5 Earth years
plt.figure(figsize= (10,8))
t = np.arange(0, 3.5, dt)
plt.plot(0,0,'*y', label = 'Sun')
plt.plot(three_body_sim(t, 'Jupiter, Earth', MJ = 1000*10E-3)[0],three_body_sim(t, 'Jupiter, Earth', MJ = 1000*10E-3)[1]         ,'--k', label = 'Jupiter')
plt.plot(three_body_sim(t, 'Jupiter, Earth', MJ = 1000*10E-3)[0][0],three_body_sim(t, 'Jupiter, Earth', MJ = 1000*10E-3)[1][0]         ,'ok')
plt.plot(three_body_sim(t, 'Jupiter, Earth', MJ = 1000*10E-3)[2],three_body_sim(t, 'Jupiter, Earth', MJ = 1000*10E-3)[3]         ,'--b', label = 'Earth')
plt.plot(three_body_sim(t, 'Jupiter, Earth', MJ = 1000*10E-3)[2][0],three_body_sim(t, 'Jupiter, Earth', MJ = 1000*10E-3)[3][0]         ,'ob')
plt.title('Orbits Over 3.5 Earth Years (MJ = 1000MJ)')
plt.xlabel('x [AU]')
plt.ylabel('y [AU]')
plt.axis('square')
plt.xlim(-10,10)
plt.ylim(-10,10)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig('Earth_Jupiter_3.5yrs.png')
plt.show()


# In[ ]:


#q2(c) Determine and plot orbits of Jupiter Asteroid system wrt. the sun over 20 Earth years
t = np.arange(0, 20, dt)
plt.figure(figsize= (10,8))
plt.plot(0,0,'*y', label = 'Sun')
plt.plot(three_body_sim(t, objects = 'Jupiter, Asteroid', MJ = 10E-3)[0],three_body_sim(t, 'Jupiter, Asteroid', MJ = 10E-3)[1]         ,'--k', label = 'Jupiter')
plt.plot(three_body_sim(t, objects = 'Jupiter, Asteroid', MJ = 10E-3)[0][0],three_body_sim(t, 'Jupiter, Asteroid', MJ = 10E-3)[1][0]         ,'ok')
plt.plot(three_body_sim(t, objects = 'Jupiter, Asteroid', MJ = 10E-3)[2],three_body_sim(t, 'Jupiter, Asteroid', MJ = 10E-3)[3]         ,'--r', label = 'Asteroid')
plt.plot(three_body_sim(t, objects = 'Jupiter, Asteroid', MJ = 10E-3)[2][0],three_body_sim(t, 'Jupiter, Asteroid', MJ = 10E-3)[3][0]         ,'or')
plt.title('Orbits Over 20 Earth Years')
plt.xlabel('x [AU]')
plt.ylabel('y [AU]')
plt.axis('square')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig('Asteroid_Jupiter_20yrs.png')
plt.show()


# In[ ]:




