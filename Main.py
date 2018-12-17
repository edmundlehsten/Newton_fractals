# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 10:02:03 2018

"""
import numpy as np

class fractal2D:
    def __init__(self, f, der = None):
        """
        f : function taking a tuple and returning a tuple
        der : reivative of function
        """
        self.func = f
        if not der is None:
            self.der = der
        else:
            der = None # TODO: need to compute derivative somehow...
        self.zeros = np.array([])
    
    def newtonMethod(self, guess):
        """
        function that carries out the newton integration method
        Author: Edmund
        Inputs
        ======
        guess - initial guess to for a zero
        
        Outputs        
        ======
        tuple - loaction of zero, returns None if the guess did not converge
        """
        maxloop = 10
        new_guess = 0 #initiate value here so that it dose not get reinitiated for each loop (to improve performance)
        for i in range(maxloop):
            new_guess = guess - self.f(guess)/self.der(guess) #should work but need a R^2 function and derivative to test...
            if np.abs(new_guess-guess) < 10**(-5): #close enough to a zero value...
                return new_guess
            guess = new_guess
        return None  # return none if did not converge
        
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
    def simpleNewtonMethod(self,guess):
        """
        Input
        =====
        guess - tuple for inital guess

        Output
        ======
        tuple, zero we converged to, None if it did not converge

        """
        maxLoop = 1000
        new_guess = np.array([0,0])
        jacob_mat = np.linalg.inv(self.partialDerivatives(guess[0],guess[1]))
        for i in range(maxLoop):
            new_guess = guess - np.matmul(jacob_mat, f(guess[0],guess[1])) # copied from task 2
            dist =  (new_guess - guess)**2
            if (dist < 10**-5).all():
                return new_guess
            guess = new_guess
        return None
