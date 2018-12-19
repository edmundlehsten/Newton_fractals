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
        self.zeros = np.array([[]])
    
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
    
    def plot(self, N, a, b, c, d):
        """
        where N eventually determines the size of the matrix and a,b,c,d are the maually
        given intervals (a,c) corresponds to the bottom left corner of the grid 
        and (b,d) the top right corner of the grid
        """
        Y, X=np.meshgrid(linspace(a,b,N),linspace(c,d,N),indexing='ij')
        """
        Gives the transpposed matrices X and Y
        """
        v_zeroes=np.vectorize(find_zero)
        """ vectorizes the function v_zeroes
        """
        A=(v_zeroes(X,Y))
        """
        creats matrix A 
        """
        pcolor(X,Y,A)
        return plt.plot(X,Y,marker='.',colour='k',linstyle='none'),
    
    """
    input: tuple
    output: matrix, matrix grid 
    """
      
     
    
    
    
    
    
    
    
    def simpleNewtonMethod(self):
        """
        (task 5)
        TODO: write method descriptiion and method...
        """
        return
    
