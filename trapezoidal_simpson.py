# Author: Milica Ivetic
# The following code are two external functions. One uses the trapezoidal 
# method of integration and the other uses Simpson's rule. 


# The following function takes code from Chapter 5.1, pg. 142 of the textbook

def trapezoidal(f,N, a, b):
    """ This function uses the Trapezoidal rule to evaluate a given function
    f from point a to b, using N slices. 
    INPUT:
    f is the function being considered in the question
    N [integer] is the number of slices 
    a, b [integer] are the starting points and end points of the integral
    
    OUTPUT:
    res [float] value of the definite integral"""
    
    h = (b-a)/N
    
    s = 0.5*f(a) +0.5*f(b)
    
    for k in range(1,N):
        s += f(a+k*h)
        
    res = h*s
    
    return res

#Using the same principles as the previous function, we write out Simpson's rule

def simpson(f, N, a, b):
    """ This function uses Simpson's rule to evaluate a given function
    f from point a to b, using N slices. 
    INPUT:
    f is the function being considered in the question
    N [integer] is the number of slices 
    a, b [integer] are the starting points and end points of the integral
    
    OUTPUT:
    res [float] value of the definite integral"""    
    
    h = (b-a)/N
    
    s = f(a) + f(b)
    
    s_odd = 0
    s_even = 0
    
    for k in range(1, N, 2):
        s_odd += f(a + k*h)
        
    for k in range(2, N, 2):
        s_even += f(a +k*h)
    
    res = (1/3)*h*(s + 4*s_odd + 2*s_even)
    
    return res