import numpy as np
import sys

from kaucherpy.core.Round import *
from kaucherpy.support.error import IntervalError
from kaucherpy.support.error import TypeIntervalError
from kaucherpy.support.error import UndefinedValueIntervalError 

__all__ = [
    "Kaucher"
]

class Kaucher(object):
    __lower = 0.0
    __upper = 0.0
    __isEmpty = False

    def __init__(self,lower=None,upper=None):
        if(upper == None and lower == None):
            self.upper = np.nan
            self.lower = np.nan
            self.isEmpty = True
        
        elif(upper is None):
            self.lower = np.float64(lower)
            self.upper = np.float64(lower)
        else:
            self.lower = np.float64(lower)
            self.upper = np.float64(upper)
        


    @property
    def lower(self):
        return self.lower
    @property
    def upper(self):
        return self.upper

    @lower.setter
    def lower(self,value):
        self.lower = np.float64(value)
        self.isEmpty = False
    
    @upper.setter
    def upper(self,value):
        self.upper = np.float64(value)
        self.isEmpty = False


    def __str__(self):
        return "["+str(self.lower)+" , "+str(self.upper)+"]"

    def __repr__(self):
        return "[%r, %r]" % (self.lower, self.upper)

    def __getitem__(self):
        return np.array([self.lower,self.upper])

    def upper(self):
        return self.upper

    def lower(self):
        return self.lower

    def setEmpty(self,v):
        self.__isEmpty = v

    def isEmpty(self,):
        return self.__isEmpty

    def __add__(self,other):
        other = self.__checkValue(other)

        Round.set_down_rounding()
        lower = self.lower+other.lower
        Round.set_up_rounding()
        Round.get_rounding()
        upper = self.upper+other.upper
        Round.set_normal_rounding()
        return Kaucher(lower,upper)



    def __sub__(self,other):
        other = self.__checkValue(other)

        Round.set_down_rounding()
        lower = self.lower-other.upper
        Round.set_up_rounding()
        upper = self.upper-other.lower
        Round.set_normal_rounding()

        return Kaucher(lower,upper)

    def __mul__(self,other):
        other = self.__checkValue(other)
        lower = np.float64(0)
        upper = np.float64(0)
        if(self.lower <= np.float64(0) and self.upper <= np.float64(0)):
            if(other.lower <= np.float64(0) and other.upper <= np.float64(0)):
                Round.set_down_rounding()
                lower = self.upper*other.upper
                Round.set_up_rounding()
                upper = self.lower*other.lower
                Round.set_normal_rounding()

            elif(other.lower < np.float64(0) < other.upper):
                Round.set_down_rounding()
                lower = self.lower*other.upper
                Round.set_up_rounding()
                upper = self.lower*other.lower
                Round.set_normal_rounding()

            elif(other.lower >= np.float64(0) and other.upper >= np.float64(0)):
                Round.set_down_rounding()
                lower = self.lower*other.upper
                Round.set_up_rounding()
                upper = self.upper*other.lower
                Round.set_normal_rounding()

            elif(other.upper < 0 < other.lower):
                Round.set_down_rounding()
                lower = self.upper*other.upper
                Round.set_up_rounding()
                upper = self.upper*other.lower
                Round.set_normal_rounding()

        elif(self.lower < np.float64(0) < self.upper):
            if(other.lower <= np.float64(0) and other.upper <= np.float64(0)):
                Round.set_down_rounding()
                lower = self.upper*other.lower
                Round.set_up_rounding()
                upper = self.lower*other.lower
                Round.set_normal_rounding()

            elif(other.lower < np.float64(0) < other.upper):
                Round.set_down_rounding()
                lower = min(self.lower*other.upper,self.upper*other.lower)
                Round.set_up_rounding()
                upper = max(self.lower*other.lower,self.upper*other.upper)
                Round.set_normal_rounding()

            elif(other.lower >= np.float64(0) and other.upper >= np.float64(0)):
                Round.set_down_rounding()
                lower = self.lower*other.upper
                Round.set_up_rounding()
                upper = self.upper*other.upper
                Round.set_normal_rounding()

            elif(other.upper < 0 < other.lower):
                lower = np.float64(0)
                upper = np.float64(0)

        elif(self.lower >= np.float64(0) and self.upper >= np.float64(0)):
            if(other.lower <= np.float64(0) and other.upper <= np.float64(0)):
                Round.set_down_rounding()
                lower = self.upper*other.lower
                Round.set_up_rounding()
                upper = self.lower*other.upper
                Round.set_normal_rounding()

            elif(other.lower < np.float64(0) < other.upper):
                Round.set_down_rounding()
                lower = self.upper*other.lower
                Round.set_up_rounding()
                upper = self.upper*other.upper
                Round.set_normal_rounding()

            elif(other.lower >= np.float64(0) and other.upper >= np.float64(0)):
                Round.set_down_rounding()
                lower = self.lower*other.lower
                Round.set_up_rounding()
                upper = self.upper*other.upper
                Round.set_normal_rounding()

            elif(other.upper < 0 < other.lower):
                Round.set_down_rounding()
                lower = self.lower*other.lower
                Round.set_up_rounding()
                upper = self.lower*other.upper
                Round.set_normal_rounding()

        elif(self.upper < 0 < self.lower):
            if(other.lower <= np.float64(0) and other.upper <= np.float64(0)):
                Round.set_down_rounding()
                lower = self.upper*other.upper
                Round.set_up_rounding()
                upper = self.lower*other.upper
                Round.set_normal_rounding()

            elif(other.lower < np.float64(0) < other.upper):
                lower = np.float64(0)
                upper = np.float64(0)
            elif(other.lower >= np.float64(0) and other.upper >= np.float64(0)):
                Round.set_down_rounding()
                lower = self.lower*other.lower
                Round.set_up_rounding()
                upper = self.upper*other.lower
                Round.set_normal_rounding()

            elif(other.upper < 0 < other.lower):
                Round.set_down_rounding()
                lower = max(self.lower*other.lower,self.upper*other.upper)
                Round.set_up_rounding()
                upper = min(self.lower*other.upper,self.upper*other.lower)
                Round.set_normal_rounding()

        return Kaucher(lower,upper)

    def __div__(self,other):
        other = self.__checkValue(other)
        if(other.lower==0 or other.upper == 0):
            raise Exception("Division by Zero")
        else:
            Round.set_down_rounding()
            lower = self.lower/other.upper
            Round.set_up_rounding()
            upper = self.upper/other.lower
            Round.set_normal_rounding()
            return Kaucher(lower,upper)



    def __or__(self, other):
        other = self.__checkValue(other)
        return Kaucher(max(self.lower,other.lower),min(self.upper,other.upper))

    def __and__(self,other):
        other = self.__checkValue(other)
        return Kaucher(min(self.lower,other.lower),max(self.upper,other.upper))

    def __invert__(self):
        return Kaucher(self.upper,self.lower)


    def __neg__(self):
        return Kaucher(-self.upper,-self.lower)

    def __pow__(self,other):
        return Kaucher(self.lower**other,self.upper**other)

    def __contains__(self,other):
        if(type(other) is not Kaucher):
            if(self.lower <= other and self.upper >= other):
                return True
            else:
                return False
        else:
            if((other.lower >= self.lower) and (self.upper>=other.upper)):
                return True
            else:
                return False

    def __eq__(self,other):
        other = self.__checkValue(other)
        if((other.lower == self.lower) and (other.upper == self.upper)):
            return True
        else:
            return False

    def __lt__(self,other):
        other = self.__checkValue(other)
        if((self.lower < other.lower) and (self.upper < other.upper)):
            return True
        else:
            return False

    def __le__(self,other):
        other = self.__checkValue(other)
        if((self.lower <= other.lower) and (self.upper <= other.upper)):
            return True
        else:
            return False

    def __gt__(self,other):
        other = self.__checkValue(other)
        if((self.lower > other.lower) and (self.upper > other.upper)):
            return True
        else:
            return False

    def __ge__(self,other):
        other = self.__checkValue(other)
        if((self.lower >= other.lower) and (self.upper >= other.upper)):
            return True
        else:
            return False


    def __checkValue(self,other):
        self.__validateDefinedValue(other)
        if(type(other) is not Kaucher):
            other = Kaucher(other)
        return other

    __radd__ = __add__
    __rsub__ = __sub__
    __rmul__ = __mul__
    __rdiv__ = __div__
    __ror__ = __or__
    __rand__ = __and__
    __rpow__ = __pow__


#a = Kaucher(2.0,3.0)
#b = Kaucher(2.5)
#d = Kaucher(5.0)
#c = 2.5
#print d > b
#print b in a
#print (c in a)
#print d in a
#print (a==b)
#print (a==a)
#print (b==2.5)
#print ~a
#print -a
#print (2.5 + b)
#print (a or b)