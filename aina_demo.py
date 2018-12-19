from scipy import *
from pylab import *
import matplotlib.pyplot as plt 
import numpy as np


def zeroes(x,y): 
 """
 pretend zero functionto test A out, has to be replaced with correct one if correct
 """
    return (sin(x**2 +y**2)/(x**2 +y**2))
def plot(N,a,b,c,d):
    """
        where N eventually determines the size of the matrix and a,b,c,d are the maually
        given intervals (a,c) corresponds to the bottom left corner of the grid 
        and (b,d) the top right corner of the grid
    """
    Y,X =np.meshgrid(linspace(a,b,N),linspace(c,d,N),indexing='ij') 
    """
    this creates the transposed matrices X and Y
    """
    v_zeroes=np.vectorize(zeroes) 
    """
    this is the vectorized function of zeroes,task 3(?)
    """
    A=(v_zeroes(X, Y) )
    """
    this creates a matrix A
    """
    pcolor(X,Y,A) 
    """probably not correct but wanted to see if A does what it should
    """
    return X, Y, plt.plot(X, Y , marker='.', color='k', linestyle='none')
plot(50,1,5,2,6)
