from scipy import *
import numpy as np
from matplotlib import pyplot as pt
from matplotlib.colors import LinearSegmentedColormap 

class fractal2D:
    def __init__(self, f, derivative = None):
        """
        Initialisation Function:
        Author: Edmund
        Inputs
        ======
        f : function taking a tuple and returning a tuple
        der : [optional] derivative of the function taking a tuple and returning a 2x2 matrix (Jacobian matrix)
        """
        self.function = f
        self.zeros = np.array([[ ]])
        if derivative==None:
            self.derivative=self.partialDerivatives
        else:
            self.derivative=derivative
    
    def partialDerivatives(self,x,y):
        """
        function that calculates the Jacobian of a DIFFERENTIABLE function
        Author: Nadine
        Inputs
        ======
        x,y - floats/ integers
        
        Outputs        
        ======
        J - Jacobian (partial derivatives) of a function (2x2 numpy array)
        """
        if not (isinstance(x,(int,float,np.int64,np.float64,np.int32,np.float32)) and isinstance(y,(int,float,np.int64,np.float64,np.int32,np.float32))):
            raise TypeError('The two input variables must be integers or floats.')
            
        f=self.function
        J=zeros(2)
        h=1e-10
        J=1/h*array([f(x+h,y)-f(x,y),f(x,y+h)-f(x,y)]).T
        return J
        
    def newtonMethod(self,guess):
        """
        function that carries out the newton integration method
        Author: Nadine
        Inputs
        ======
        guess - (tuple), list or 1-dimensional array of length 2 with float/integer entries
        
        Outputs        
        ======
        zero - loaction of zero, returns None if the guess did not converge (tuple)
        """
        #check if inputs are correct
        if isinstance(guess,(list,tuple,np.ndarray)) and size(guess)==2:
            guess=array(guess)
        else:
            if not isinstance(guess,(list,tuple,np.ndarray)):
                raise TypeError('The initial guess is not of the correct type.\nIt is of type:' + str(type(guess)))
            else     
                raise TypeError('The initial guess is not of the correct length.\nIt is of length:' + str(size(guess)))

        #run simple newton method
        maxloop = 1000
        new_guess = array([0,0]) #initiate to reduce constant reinitiation
        for i in range(maxloop):
            Jacobian=self.derivative(guess[0],guess[1])
            if np.linalg.det(Jacobian)==0:
                return None,maxloop
            else:
                Jacobian_inv=np.linalg.inv(Jacobian)
            new_guess = guess - np.matmul(Jacobian_inv,self.function(guess[0],guess[1]))

            #test if as a stationary point -> at zero
            if np.abs(new_guess[0]-guess[0]) < 10**(-5) and np.abs(new_guess[1]-guess[1]) < 10**(-5): 
                return new_guess , i
            guess = new_guess
        else:
            return None,maxloop # return None if did not converge
 

    def find_zero(self, guess, simple_newton = 0):
        """
        function that carries out the newton integration method
        Author: Edmund
        Inputs
        ======
        guess - initial guess to for a zero
        simple_newton - [optional] used to define which newton method to use
            0 - normal method
            1 - simple method 
            2 - optimized method (combination of 0 & 1)
        
        Output
        ======
        returns index of the found zero, None if point did not converge
        """
        val , i = self.newton(guess)
        tol = 10**(-4)  # tolerance value...
        if val is None: # value did not converge
            return -1,i
        if self.zeros.size == 0: # first zero using the current method 
            self.zeros = np.array([val])

        t = (self.zeros-val)**2
        dist = t[:,0] + t[:,1]
        if (dist<tol).any() : #zero already exists
            return np.where(dist<tol)[0][0] , i
        else:  # value dose not exist yet
            self.zeros=np.reshape(np.append(self.zeros,val),(-1,2)) #add value to zeros
            return self.zeros.size/2-1 , i 
    
    def call_findZero(self,A,B,simple = False):
        """
        a helper function such we can vectorise without an error,
        function takes two inputs and then calls find zero combining hte two inputs into one tuple
        """
        return (self.find_zero((A,B),simple)) 

    def plot(self, N, a, b, c, d, simple_newton = 0, show_plot = True):
        """
        Function used to plot fractal within a certain range:
        Author: Aina
        Inputs:
        ======
        N - (int) number of pixels in each direction
        a - (int) lower x bound
        b - (int) upper x bound
        c - (int) lower y bound 
        d - (int) upper y bound
        simple_newton - [optional] (int) which newton method to use:
            0 - normal newton approximation
            1 - simple newton approximation
            2 - optimised approcimation using both (0 & 1)
        show_plot - [optional] (boolean) wether or not to show the generated plot

        outputs:
        ========
        A - a NxN matrix containing int indecies one for each zero, -1 if the point did not converge
        B - a index containing how many steps it took for the newton method to determine convergence
        """
        #implementing simplified newton, implemented in this way to speed up computations
        if simple_newton == 1:                                                  
                 self.newton =  self.simpleNewtonMethod                              
        elif simple_newton == 2:                                                
             self.newton = self.optimised
        else:
            self.newton =  self.newtonMethod

        # Generate a matrix containing X coordinates and one containing Y coordinates
        # for each point in the NxN matrix to return
        Y, X=np.meshgrid(linspace(a,b,N),linspace(c,d,N),indexing='ij')

        # Vectorising the function find zero to run it on the entire matrix
        v_zeroes=np.vectorize(self.call_findZero)

        # Generate A,B matrix
        A, B=(v_zeroes(X,Y,simple_newton))
        # if plot is to be shown (following code by Edmund):
        if show_plot:
            # Addjust B values to make for nicer graph
            B = ((B-1)%10)/10
            # get set1 color map, bright easily distingueshed colors:
            color_map = pt.cm.get_cmap('Set1')
            # set color for values less than 0 to black
            color_map.set_under('0')
            # plot the A matrix (indecies)
            pcol = pt.pcolormesh(X,Y,A,cmap=pt.cm.Set1,edgecolor = 'face', linewidth=0,rasterized=True, vmin=0,vmax = A.max())
            # display color bar at the side
            pt.colorbar() 
            #make custom color map black backgound changing alpha gradient
            custom = [(0,(0, 0, 0, 0)), (1,(0, 0, 0, 0.5))]
            my_cmap = LinearSegmentedColormap.from_list('something',custom)
            # plot the gradient colormap generated on top of colored plot
            pcol = pt.pcolormesh(X,Y,B,edgecolor = 'face',cmap=my_cmap,linewidth=0,rasterized=True)
            #show the plot
            pt.show()
        return A,B
    

    def simpleNewtonMethod(self,guess):
        """
        function that carries out the simple newton integration method
        Author: Edmund
        Inputs
        ======
        guess - tuple, list or 1-dimensional array of length 2 with float/ integer entries
        
        Outputs        
        ======
        zero - loaction of zero, returns None if the guess did not converge (tuple)
        """
        # check that input guess value is of the correct type
        if isinstance(guess,(list,tuple,np.ndarray)) and size(guess)==2:
            guess=np.array(guess)
        else:
            raise TypeError('The initial guess is not of the correct type. It must be a list with two integer/ float entries, a tuple or an array of size (2,)')
        # run simple method
        maxloop = 1000
        new_guess = np.array([0,0]) #initiate value to not do it repeatedly in loop 
        # Compute jacobian matrix
        Jacobian=self.derivative(guess[0],guess[1])
        #check that it is invertible and invert
        if np.linalg.det(Jacobian)==0:  
            return None,maxloop
        else:
            Jacobian_inv=np.linalg.inv(Jacobian)
        # itterativly approach zero point 
        for i in range(maxloop):
            new_guess = guess - np.matmul(Jacobian_inv,self.function(guess[0],guess[1]))
            if (np.abs(new_guess)>1e+30).any(): #check if we are running into an overflow errer => diverges
                return None,maxloop
            if np.abs(new_guess[0]-guess[0]) < 10**(-5) and np.abs(new_guess[1]-guess[1]) < 10**(-5): #converged!
                return new_guess , i
            guess = new_guess
        else:
            return None,maxloop # return None if did not converge
        
    def optimised(self,guess):                                                  
        """
        function that uses both simple and normal newton method to try and get the best of both.
        Author: Edmund
        Inputs:
        ======
        guess - (tuple, len2) initial starting point

        Output:
        ======
        returns simple newton zero and iterations if it converged, otherways returns the normal newton zero and number of itterations
        """
        ret_val, i = self.simpleNewtonMethod(guess)   
        if isinstance(ret_val,(np.ndarray)):                                                    
            return ret_val, i                                                   
        else :                                                                  
            return self.newtonMethod(guess) 
            
#TODO:------------------THIS ALL NEEDS TO GO!!!---------------------------------
def f(x,y):
    return np.array([x**3-3*x*y**2-1,3*x**2*y-y**3])
k = fractal2D(f)
#k.plot(250,-0.5,0.5,-0.5,0.5,True)
def diff(guess):
    print("Diff:", k.simpleNewtonMethod(guess)-k.newtonMethod(guess))

def g(x,y):
    return np.array([(x**3)-(3*x)*(y**2)-(2*x)-2,3*(x**2)*y-(y**3)-2*y])
l=fractal2D(g)

def h(x,y):
    return np.array([x**8-28*x**6*y**2+70*x**4*y**4+15*x**4-28*x**2*y**6-90*x**2*y**2+y**8+15*y**4-16,8*x**7*y-56*x**5*y**3+56*x**3*y**5+60*x**3*y-8*x*y**7-60*x*y**3])
m=fractal2D(h)
