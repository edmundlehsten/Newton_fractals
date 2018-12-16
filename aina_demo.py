from scipy import *
from pylab import *
import matplotlib.pyplot as plt 
import numpy as np


def zeroes(x,y): ##pretend zero functionto test A out, has to be replaced with correct one if correct
    return (sin(x**2 +y**2)/(x**2 +y**2))
def plot(N,a,b,c,d):
    #function plot has N,a,b,c,d as parameters. N eventually deterines the shape of the matrices.
    #a,b,c,d are manually given numbers as intervals
    X,Y =np.meshgrid(linspace(a,b,N),linspace(c,d,N),indexing='ij') 
    #this creates the matrices X and Y
    X=X.transpose() #This transposes them
    Y=Y.transpose()
    v_zeroes=np.vectorize(zeroes) #this is the vectorized function of zeroes,task 3(?)
    A=(v_zeroes(X, Y) )#this creates a matrix A
    print(pcolor(X,Y,A)) #probably not correct but wanted to see if A does what it should
    return X, Y, plt.plot(X, Y , marker='.', color='k', linestyle='none')
plot(50,1,5,2,6)
