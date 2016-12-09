from astromodels import Function1D, FunctionMeta
import astropy.units as u
import astropy.constants as constants
from scipy.integrate import quad

import numpy as np



# Define astropy models

class KRL(Function1D):
    r"""
    description :
        The Pulse function from Kocevski, Ryde, and Liang (2003)
    latex : $  $
    parameters :
        K :
            desc : normalization
            initial value : 1
            min : 0
    
        t_max :
            desc : time of maximum flux
            initial value : 5
            min : 0.
    
  
        rise:
            desc : time rise constant
            initial value : 1.
            min : 0
            
            
            
        decay:
            desc : time decay constant
            initial value : 2.5
            min : 0
            fix: no
            
        c:
            desc : shift factor
            initial value : 0
            fix: yes
       
    """

    __metaclass__ = FunctionMeta

    def _set_units(self, x_unit, y_unit):

       
        self.K.unit = y_unit
        self.t_max.unit = x_unit
        self.c.unit = x_unit
        self.rise.unit = u.dimensionless_unscaled
        self.decay.unit = u.dimensionless_unscaled
      


    def evaluate(self, x, K, t_max, rise, decay, c):


        f = K*(np.power((((x+c)/(t_max+c))),rise)/power(((d+(rise*np.power((((x+c)/(t_max+c))),(1+rise))))/(decay+r)),((decay+rise)/(1+rise))))

        return f
        



class Norris(Function1D):
    r"""
    description :
        The Pulse function from Kocevski, Ryde, and Liang (2003)
    latex : $  $
    parameters :
        K :
            desc : normalization
            initial value : 1
            min : 0
    
        t_start :
            desc : start time 
            initial value : 5
            min : 0.
    
  
        t_rise:
            desc : time rise constant
            initial value : 1.
            min : 0
            
            
            
        t_decay:
            desc : time decay constant
            initial value : 2.5
            min : 0
            fix: no
            
        
       
    """

    __metaclass__ = FunctionMeta

    def _set_units(self, x_unit, y_unit):

       
        self.K.unit = y_unit
        self.t_start.unit = x_unit
        self.t_rise.unit = x_unit
        self.t_decay.unit = x_unit
        
      


    def evaluate(self, x, K, t_start, t_rise, t_decay):



       f = K*np.exp(2*(t_rise/ t_decay)**(1/2) ) * np.exp( -t_rise / (x - t_start) - (x - t_start) / t_decay )


       return f
        
