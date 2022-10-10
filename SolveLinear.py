# SolveLinear.py
# Python module for PHY407
# Paul Kushner, 2015-09-26
# Modifications by Nicolas Grisouard, 2018-09-26
# This module contains useful routines for solving linear systems of equations.
# Based on gausselim.py from Newman
#from numpy import empty
# The following will be useful for partial pivoting
from numpy import empty, copy
import numpy as np

# PartialPivot MODIFIED BY: MILICA IVETIC 2022-10-04

def GaussElim(A_in, v_in):
    """Implement Gaussian Elimination. This should be non-destructive for input
    arrays, so we will copy A and v to
    temporary variables
    IN:
    A_in, the matrix to pivot and triangularize
    v_in, the RHS vector
    OUT:
    x, the vector solution of A_in x = v_in """
    # copy A and v to temporary variables using copy command
    A = copy(A_in)
    v = copy(v_in)
    N = len(v)

    for m in range(N):
        # Divide by the diagonal element
        div = A[m, m]
        A[m, :] /= div
        v[m] /= div

        # Now subtract from the lower rows
        for i in range(m+1, N):
            mult = A[i, m]
            A[i, :] -= mult*A[m, :]
            v[i] -= mult*v[m]

    # Backsubstitution
    # create an array of the same type as the input array
    x = empty(N, dtype=v.dtype)
    for m in range(N-1, -1, -1):
        x[m] = v[m]
        for i in range(m+1, N):
            x[m] -= A[m, i]*x[i]
    return x


def PartialPivot(A_in, v_in):
    """ In this function, code the partial pivot (see Newman p. 222) """
    A = copy(A_in)
    v = copy(v_in)
    N = len(v)

    for m in range(N):
        
        if np.abs(A[m,m]) < 1.0e-16: #if the element in matrix is ~0
            for i in range(m+1, N): #considering lower rows
                if np.abs(A[i,m]) > np.abs(A[m,m]):
                    #swapping rows and corresponding vector elements
                    A[m,:], A[i,:] = copy(A[i,:]), copy(A[m,:])
                    v[m], v[i] = copy(v[i]), copy(v[m])
        
        
        #The rest is same as GaussElim
        # Divide by the diagonal element
        div = A[m, m]
        A[m, :] /= div
        v[m] /= div

        # Now subtract from the lower rows
        for i in range(m+1, N):
            mult = A[i, m]
            A[i, :] -= mult*A[m, :]
            v[i] -= mult*v[m]

    # Backsubstitution
    # create an array of the same type as the input array
    x = empty(N, dtype=v.dtype)
    for m in range(N-1, -1, -1):
        x[m] = v[m]
        for i in range(m+1, N):
            x[m] -= A[m, i]*x[i]
    return x    
    

#IGNORE

#Atest2 = np.array([[2, 1, 4, 1], [3,4,-1, -1],[1,-4,1, 5],[2,-2,1, 3]], float)
#Vtest = np.array([1, 0], float)

#Atest = np.array([[1e-20, 1], [1, 1]], float)
#Vtest2 = np.array([-4, 3, 9,7], float)

#print(GaussElim(Atest, Vtest))
#print(PartialPivot(Atest, Vtest))

#print(GaussElim(Atest2, Vtest2))
#print(PartialPivot(Atest2, Vtest2))