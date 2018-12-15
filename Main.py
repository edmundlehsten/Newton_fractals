import numpy as np
class fractal2D:
    def __init__(self, f, derivative = None):
        """
        f : function taking a tuple and returning a tuple
        der : derivative of the function taking a tuple and returning a 2x2 matrix (Jacobian matrix)
        """
        self.function = f
        self.zeros = np.array([])
        if derivative==None:
            self.derivativeInput=str(None)
        else:
            self.derivativeInput=derivative
       
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
        if not (isinstance(x,(int,float)) and isinstance(y,(int,float))):
            raise TypeError('The two input variables must be integers or floats.')
            
        if self.derivativeInput=="None":
            f=self.function
            J=zeros(2)
            h=1e-10
            J=1/h*array([f(x+h,y)-f(x,y),f(x,y+h)-f(x,y)]).T
            return J
        else:
            return self.derivativeInput(x,y)
    
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
        if isinstance(guess,(list,tuple)) and size(guess)==2:
            guess=array(guess)
        elif not (str(type(guess))=="<class 'numpy.ndarray'>" and size(guess)==2):
            raise TypeError('The initial guess is not of the correct type. It must be a list with two integer/ float entries, a tuple or an array of size (2,)')
        
        maxloop = 1000
        new_guess = array([0,0]) #initiate value here so that it dose not get reinitiated for each loop (to improve performance)
        for i in range(maxloop):
            new_guess = guess - np.matmul(np.linalg.inv(self.partialDerivatives(guess[0],guess[1])),f(guess[0],guess[1]))#should work but need a R^2 function and derivative to test...
            if np.abs(new_guess[0]-guess[0]) < 10**(-5) and np.abs(new_guess[1]-guess[1]) < 10**(-5): #close enough to a zero value...
                return new_guess
            guess = new_guess
        else:
            return None  # return None if did not converge
        

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
    