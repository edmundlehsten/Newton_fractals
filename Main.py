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
        #print(type(x))  #if there is an error message here, print x to check the type
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
            new_guess = guess - np.matmul(np.linalg.inv(self.derivative(guess[0],guess[1])),self.function(guess[0],guess[1]))#should work but need a R^2 function and derivative to test...

            if np.abs(new_guess[0]-guess[0]) < 10**(-5) and np.abs(new_guess[1]-guess[1]) < 10**(-5): #close enough to a zero value...
                return new_guess
            guess = new_guess
        else:
            return None  # return None if did not converge

  
    
  
    def find_zero(self, guess):
        """
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

