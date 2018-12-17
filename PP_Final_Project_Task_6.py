# -*- coding: utf-8 -*-
"""
Created on Sun Dec 16 21:13:55 2018

@author: Nadine
"""
from scipy import *
import numpy as np

class fractal2D:
    def partialDerivatives(self,x,y):
        """
        function that calculates the Jacobian of a DIFFERENTIABLE function
        Author: Nadine
        Inputs
        ======
        x,y - floats/ integers
        
        Outputs        
        ======
        J - Jacobian (partial derivatives) of a function
        """
        #print(type(x))
        if not (isinstance(x,(int,float,np.int64,np.float64,np.int32,np.float32)) and isinstance(y,(int,float,np.int64,np.float64,np.int32,np.float32))):
            raise TypeError('The two input variables must be integers or floats.')
            
        f=self.function
        J=zeros(2)
        h=1e-10
        J=1/h*array([f(x+h,y)-f(x,y),f(x,y+h)-f(x,y)]).T
        return J
        
    def __init__(self, f, derivative = None):
        """
        f : function taking a tuple and returning a tuple
        der : derivative of the function taking a tuple and returning a 2x2 matrix (Jacobian matrix)
        """
        self.function = f
        self.zeros = np.array([])
        if derivative==None:
            self.derivative=self.partialDerivatives
        else:
            self.derivative=derivative
    
    def newtonMethod(self,guess):
        """
        function that carries out the newton integration method
        Author: Edmund
        Inputs
        ======
        guess - tuple, list or 1-dimensional array of length 2 with float/ integer entries
        
        Outputs        
        ======
        zero - loaction of zero, returns None if the guess did not converge (tuple)
        """
        if isinstance(guess,(list,tuple,np.ndarray)) and size(guess)==2:
            guess=array(guess)
        else:
            raise TypeError('The initial guess is not of the correct type. It must be a list with two integer/ float entries, a tuple or an array of size (2,)')
        
        maxloop = 1000
        new_guess = array([0,0]) #initiate value here so that it dose not get reinitiated for each loop (to improve performance)
        for i in range(maxloop):
            new_guess = guess - np.matmul(np.linalg.inv(self.derivative(guess[0],guess[1])),f(guess[0],guess[1]))#should work but need a R^2 function and derivative to test...
            if np.abs(new_guess[0]-guess[0]) < 10**(-5) and np.abs(new_guess[1]-guess[1]) < 10**(-5): #close enough to a zero value...
                return new_guess
            guess = new_guess
        else:
            return None  # return None if did not converge
        

def f(x,y):
    return array([x**3-3*x*y**2-1,3*x**2*y-y**3])

def fp(x,y):
    return array([[3*x**2-3*y**2,-6*x*y],[6*x*y,3*x**2-3*y**2]])

new=fractal2D(f)
new1=fractal2D(f,fp)