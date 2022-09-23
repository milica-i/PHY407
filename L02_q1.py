#!/usr/bin/env python
# coding: utf-8

# In[1]:


""" Written by Madeline Nardin September 2022 """
import numpy as np
import matplotlib.pyplot as plt


# In[2]:


#Q1(a) Pseudocode 
#i. Read in supplied data set 'cdata.txt'
#ii. Define Method 1 function for eq. 1 sigma = sqrt((1/(n-1))Σ(x_i-x_bar)^2) and compute with supplied data set
#iii. Define Method 2 function for eq. 2 sigma = sqrt((1/(n-1))Σx_i^2-n*x_bar^2) and compute with supplied data set
#iv. Determine correct value for std by numpy.std()
#v. Determine relative error of eq1 & eq2 wrt. numpy.std() function and determine the more accurate method 
    #(smallest error)


# In[29]:


#load supplied data set
data = np.loadtxt('cdata.txt')


# In[12]:


#Define functions 
def eq1(dataset):
    """ Determines the standard deviation with a two pass method.

            Determines the standard deviation with a two pass method with the following equation\
            sigma = sqrt((1/(n-1))Σ(x_i-x_bar)^2).
        Parameters
        --------
        dataset : array
            dataset is any given array of numbers.

        Returns
        -------
        sigma : float
            Standard deviation of the dataset.
    """
    n = len(dataset)
    x_bar = (1/n)*np.sum(dataset)
    for i in dataset:
        std_sqd = (1/(n-1))*np.sum((i-x_bar)**2)
    if std_sqd < 0:
        return 'Negative squareroot object - imaginary standard deviation'
    else: 
        sigma = np.sqrt(std_sqd)
        print('Method 1 solution = ', sigma)
        return sigma

def eq2(dataset):
    """ Determines the standard deviation with a single pass method.

            Determines the standard deviation with a single pass method with the following equation\
            sigma = sqrt((1/(n-1))Σx_i^2-n*x_bar^2).
        Parameters
        --------
        dataset : array
            dataset is any given array of numbers.

        Returns
        -------
        sigma : float
            Standard deviation of the dataset.
    """
    n = len(dataset)
    for i in dataset:
        std_sqd = (1/(n-1))*np.sum(i**2-n*(i/n)**2)
    if std_sqd < 0:
        return 'Negative squareroot object - imaginary standard deviation'
    else: 
        sigma = np.sqrt(std_sqd)
        print('Method 2 solution = ', sigma)
        return sigma
    
def relative_error(x,x_true):
    """ Determines the relative error between two values.

            Determines the relative error between two values with the following equation\
            error = ((x-x_true)/x_true).
        Parameters
        --------
        x : float
            x is the value in question.
        
        x_true : float
            x_true is the value considered to be correct.

        Returns
        -------
        err : float
            Relative error of value in question, x, with respect to the correct value x_true.
    """
    err = ((x-x_true)/x_true)
    print('Relative err = ', err)
    return err


# In[13]:


#Compute standard deviation of provided dataset with method 1, method 2 and real value (numpy.std())
sol_eq1 = eq1(data)
sol_eq2 = eq2(data)
sol = np.std(data, ddof = 1)
print('Real solution = ', sol)


# In[14]:


#Determine relative error for method 1 and method 2 wrt. real value (numpy.std())
sol_eq1_err = relative_error(float(sol_eq1),sol)
sol_eq2_err = relative_error(float(sol_eq2),sol)

#Determine which method is more accurate (smaller error)
if sol_eq1_err < sol_eq2_err: print('Method 1 is more accurate')
else:
    print('Method 1 is more accurate')


# In[15]:


#Define sequence a and b as described in writeup 
a = np.random.normal(0.0, 1.0, 2000)
b = np.random.normal(1.0e7, 1.0, 2000)


# In[22]:


#Compute standard deviation of sequence a with method 1, method 2 and real value (numpy.std())
a_eq1 = eq1(a)
a_eq2 = eq2(a)
a_sol = np.std(a, ddof = 1)
print('Real solution = ', a_sol)


# In[23]:


#Compute standard deviation of sequence a with method 1, method 2 and real value (numpy.std())
b_eq1 = eq1(b)
b_eq2 = eq2(b)
b_sol = np.std(b, ddof = 1)
print('Real solution = ', b_sol)


# In[36]:


#Determine relative error for sequence a method 1 and method 2 wrt. real value (numpy.std())
a_sol1_err = relative_error(float(a_eq1),a_sol)
a_sol2_err = relative_error(float(a_eq2),a_sol)

#Determine which method is more accurate (smaller error) for sequence a
if a_sol1_err < a_sol2_err: print('Method 1 is more accurate for sigma(a)')
else:
    print('Method 1 is more accurate for simga(a)')


# In[37]:


#Determine relative error for sequence b method 1 and method 2 wrt. real value (numpy.std())
b_sol1_err = relative_error(float(b_eq1),b_sol)
b_sol2_err = relative_error(float(b_eq2),b_sol)

#Determine which method is more accurate (smaller error) for sequence b
if b_sol1_err < b_sol2_err: print('Method 1 is more accurate for sigma(b)')
else:
    print('Method 1 is more accurate for sigma(b)')


# In[73]:


#note that the sequence with smaller mean optimized performace of single pass method 

#Q1(d) Pseudocode 
#i. Define Method 2 function for eq. 2 sigma = sqrt((1/(n-1))Σx_i^2-n*x_bar^2) 
    #with i shifted negatively by  first index of the data set - this is 
    #done to try to shift the mean to zero while maintaining the original
    #determine the more accurate method standard deviation.
#ii. Compute with sequence a and b
#iii. Determine relative error of revised eq2 wrt. numpy.std() function and
    #determine (smallest error) 
#iv. Compare results to original method 2.


# In[71]:


def eq2_revised(dataset):
    """ Determines the standard deviation with a single pass method.

            Determines the standard deviation with a single pass method with the following equation\
            sigma = sqrt((1/(n-1))Σx_i^2-n*x_bar^2).
        Parameters
        --------
        dataset : array
            dataset is any given array of numbers.

        Returns
        -------
        sigma : float
            Standard deviation of the dataset.
    """
    n = len(dataset)
    for i in dataset:
        x = i-dataset[0]
        std_sqd = (1/(n-1))*np.sum(x**2-n*(x/n)**2)
    if std_sqd < 0:
        return 'Negative squareroot object - imaginary standard deviation'
    else: 
        sigma = np.sqrt(std_sqd)
        print('Method 2 solution = ', sigma)
        return sigma


# In[65]:


a_eq2_rev = eq2_revised(a)
a_sol2_rev_err = relative_error(float(a_eq2_rev),b_sol)


# In[66]:


#Determine which method is more accurate (smaller error)
if a_sol2_rev_err < a_sol2_err: print('Revised Method 2 is more accurate than Method 2')
if a_sol2_rev_err < a_sol1_err : print('Revised Method 2 is more accurate than Method 1')
else: print('Method 1 is more accurate than Revised Method 2')


# In[72]:


b_eq2_rev = eq2_revised(b)
b_sol2_rev_err = relative_error(float(b_eq2_rev),b_sol)


# In[68]:


#Determine which method is more accurate (smaller error)
if b_sol2_rev_err < b_sol2_err: print('Revised Method 2 is more accurate than Method 2')
if b_sol2_rev_err < b_sol1_err : print('Revised Method 2 is more accurate than Method 1')
else: print('Method 1 is more accurate than Revised Method 2')

