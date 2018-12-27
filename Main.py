from scipy import *
import numpy as np
from matplotlib import pyplot as pt

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
        self.zeros = np.array([[ ]])
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
                return new_guess , i
            guess = new_guess
        else:
            return None  # return None if did not converge

  
    
  
    def find_zero(self, guess, simple_newton = False):
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
        val , i = self.newton(guess)
        tol = 10**(-5)  # tolerance value...
        if val is None: # value did not converge
            return-1 
        if self.zeros.size == 0:
            self.zeros = np.array([val])
        t = (self.zeros-val)**2
        dist = t[:,0] + t[:,1]
        if (dist<tol).any() : #zero exists
            return np.where(dist<tol)[0][0] , i
        else:  # value dose not exist yet
            np.reshape(np.append(self.zeros,val),(-1,2)) #add value to zeros
            return self.zeros.size-1 , i 
    
    def call_findZero(self,A,B,simple = False):
        """
        a helper function such we can vectorise without an error
        """
        #self.curr += 1
        #print(self.curr, " out of ", self.max)
        return (self.find_zero((A,B),simple)) 

    def plot(self, N, a, b, c, d, simple_newton = False):
        """
        where N eventually determines the size of the matrix and a,b,c,d are the maually
        given intervals (a,c) corresponds to the bottom left corner of the grid 
        and (b,d) the top right corner of the grid
        """
        #implementing simplified newton, implemented in this way to speed up computations
        if simple_newton == True:
            self.newton =  self.simpleNewtonMethod
        else:
            self.newton =  self.newtonMethod
        #used for tracking progress...
        #self.max = N**2
        #self.curr = 0

        Y, X=np.meshgrid(linspace(a,b,N),linspace(c,d,N),indexing='ij')
        """
        Gives the transpposed matrices X and Y
        """
        v_zeroes=np.vectorize(self.call_findZero)
        """ vectorizes the function v_zeroes
        """
        A, B=(v_zeroes(X,Y,simple_newton))
        """
        creats matrix A 
        """
        pcolor(X,Y,A)
        #return plt.plot(X,Y,marker='.',colour='k',linstyle='none'),
        return A
    

    def simpleNewtonMethod(self,guess):
        """
        Input
        =====
        guess - tuple for inital guess

        Output
        ======
        tuple, zero we converged to, None if it did not converge
        """

        initial_guess = guess
        maxLoop = 100000 #need a lot more to get decent results :(
        new_guess = np.array([0,0])
        jacob_mat = np.linalg.inv(self.partialDerivatives(guess[0],guess[1])) / 2
        f = self.function(guess[0],guess[1])
        pre_dist = None
        for i in range(maxLoop):
            new_guess = guess - np.matmul(jacob_mat,self.function(guess[0],guess[1])) # copied from task 2
            print(new_guess)
            dist =(new_guess - guess)**2
            if (dist == float("inf")).any :
                return None
            if dist[0]+dist[1]<10**-8 :
                return new_guess
            guess = new_guess
            if pre_dist == None:
                pre_dist = dist[0]+dist[1]
            elif pre_dist <  dist[0]+dist[1]:
                gone_up = True
            pre_dist = dist[0]+dist[1]

        if gone_up == False:
            return new_guess
        return None

def f(x,y):
    return np.array([x**3-3*x*y**2-1,3*x**2*y-y**3])
k = fractal2D(f)
#k.plot(250,-0.5,0.5,-0.5,0.5,True)
def diff(guess):
    print("Diff:", k.simpleNewtonMethod(guess)-k.newtonMethod(guess))
