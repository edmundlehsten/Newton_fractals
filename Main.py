# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 10:02:03 2018

"""
import numpy as np

class fractal2D:
    def __init__(self, f, derivative = None):
        """
        f : function taking a tuple and returning a tuple
        der : derivative of the function taking a tuple and returning a 2x2 matrix (Jacobian matrix)
        """
        self.function = f
        if not derivative is None:
            self.derivative = derivative
        else:
            derivative = None # TODO: need to compute derivative somehow...
        self.zeros = np.array([])
    
    def newtonMethod(self,guess):
        """
        function that carries out the newton integration method
        Author: Edmund
        Inputs
        ======
        guess - tuple
        
        Outputs        
        ======
        zero - loaction of zero, returns None if the guess did not converge (tuple)
        """
        #if (not (isinstance(guess,array) and (size(guess)==2)) or (not (isinstance(guess,list) and size(guess)==2 and isinstance(guess[0],(int,float)) and isinstance(guess,(int,float))))):
        #    raise TypeError('The given initial value is not of the right type')
            
        maxloop = 1000
        new_guess = array([0,0]) #initiate value here so that it dose not get reinitiated for each loop (to improve performance)
        for i in range(maxloop):
            new_guess = guess - np.matmul(np.linalg.inv(self.derivative(guess[0],guess[1])),self.function(guess[0],guess[1])) #should work but need a R^2 function and derivative to test...
            if np.abs(new_guess[0]-guess[0]) < 10**(-5) and np.abs(new_guess[1]-guess[1]) < 10**(-5): #close enough to a zero value...
                return new_guess
            guess = new_guess
        else:
            return None  # return none if did not converge

# %% new cell
        
    def find_zero(self, guess):
        """
        function that carries out the newton integration method
        Author: Edmund
        Inputs
        ======
        guess - initial guess to for a zero
        
        Output
        ======
        returns index of zero, None if point did not converge
        """
        val = newtonMethod(guess)
        if val is None: # value did not converge
            return None        
        if (self.zeros-val < 10**-5).any: #zero already exists
            np.where(self.zeros-val < 10**-5)
        else:  # value dose not exist yet
            np.append(self.zeros,val) #add value to zeros
            return self.zeros.size-1
    
    def plot(self):
        """
        (Task 4)
        TODO: Write description of method... 
        """
        return
    def simpleNewtonMethod(self):
        """
        (task 5)
        TODO: write method descriptiion and method...
        """
        return
    