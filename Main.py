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
    def simpleNewtonMethod(self,N,a,b,c,d):
        #gen grid matrix
        X,Y = np.meshgrid(np.linspace(a,c,N),np.linspace(b,d,N))
        # def method for one step at point x,y
        def approximate(x,y):
            return np.array([x,y]) - (self.f(x,y)*self.der(x,y).inverse())  
            # not tested, probably wront sintex but the general calculation should be correct
        #vectorise previous method
        approx = np.vectorise(approximate)
        #carry out vectorised method on entire grid
        approx_mat = approx(X,Y)
        # def get index
        max_index = 0
        def get_index(x,y):
            if index[x,y]>-1: # if we already know the value here
                return index[x,y]
            delta = approx_mat[x,y]
            new_x = delta[0]
            new_y = delta[1]
            if new_x > np.max(a,c) or new_x < np.min(a,c):
                return 0  # out of bounds for x hence we gues it dose not converge...
            elif new_y < np.min(b,d) or new_y > np.max(a,b):
                return 0  # out of bounds for y hence we guess it dose not converge...
            
            mapx = int((N*(new_x-a)/(c-a))+0.5)  # add 0.5 to improve rounding and prevent it from always rounding down.
            mapy = int((N*(new_y-b)/(d-b))+0.5)

            if mapx == x and mapy == y:
                #on same spot -> at zero -> return next index
                max_index += 1
                index[x,y] = max_index
            index[x,y] = get_index(mapx,mapy)
            return index[x,y] 

        #vectorise previous definition
        get_index_mat = np.vectorise(get_index) 
        #gen index matrix to store where each point converged to
        index = -1 * np.ones(N,N,dtype = int)
        #run vectorised matrix on whole index matrix
        index = get_index_mat(X,Y)
        #return index matrix
        return index
