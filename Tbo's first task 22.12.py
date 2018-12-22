# -*- coding: utf-8 -*-
"""
Created on Sat Dec 22 10:32:12 2018

@author: thiba
"""

class ErrorLenFunction(Exception) :
    pass
      
class ErrorDimFunction(Exception) :
    pass
   
class ErrorTypeOutFunct(Exception) :
    pass

class ErrorLenDerivative(Exception) :
    pass
  
class ErrorDimDerivative(Exception) :
    pass
             
class ErrorTypeDerivative(Exception) :
    pass



def f(x,y):
    if  len(str(signature(f)).split(",")) != 2 :
        raise ErrorLenFunction("f has to take 2 arguments")
    return   3*x + 7**y, 3*x + 2*y



class fractal2D:
    def __init__(self, f, derivative = None):
        """
        f : function taking a tuple and returning a tuple
        der : derivative of the function taking a tuple and returning a 2x2 matrix (Jacobian matrix)
        """
        
        """
        Check The {len/type/dimension} of the {Input/Output} of {f and der}
        """
        self.function = f
        self.zeros = np.array([])
            
        if  len(str(signature(f)).split(",")) != 2 :
            raise ErrorLenFunction("f has to take 2 arguments")
        returned = f(1,1)
        if not isinstance(returned,(list,tuple,np.ndarray)):
            raise ErrorTypeOutFunct("f is of wrong return type")
        if len(returned) != 2 :
            raise ErrorDimFunction("f returned wrong dimension")
        if not isinstance(returned[0],(float,int,np.float64, np.float32, np.int32, np.int64)) or not isinstance(returned[1],(float,int,np.float64, np.float32, np.int32, np.int64)) :
            raise ErrorTypeOutFunct("f returned wrong value type")
        
        
        if derivative != None :            
            if len(str(signature(derivative)).split(",")) != 2 :
                raise ErrorLenDerivative("f has to take 2 arguments")
            if not isinstance(derivative,(list,tuple,np.ndarray)):
                raise ErrorTypeDerivative("derivative is of wrong return type")
            if len(derivative) != 2 :
                raise ErrorDimDerivative("derivative returned wrong dimension")
            if not isinstance(derivative[0],(float,int,np.float64, np.float32, np.int32, np.int64)) or not isinstance(derivative[1],(float,int,np.float64, np.float32, np.int32, np.int64)) :
                raise ErrorTypeDerivative("derivative returned wrong value type")
