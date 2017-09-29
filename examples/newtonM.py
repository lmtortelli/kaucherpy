''' Newton_method3.py
find value for x in f(x) = 0 using Newton's method
see:
http://en.wikipedia.org/wiki/Newton%27s_method
tested with Python27, IronPython273 and Python33 by vegaseat 08jan2013
'''
from kaucherpy import *
from kaucherpy.core import *
from kaucherpy.kaucher import *
from intpy import *

def powI(xInt,x):
    return IReal(xInt.inf**x,xInt.sup**x)

def derivative(f):
    def compute(x, dx):
        return (f(x+dx) - f(x))/dx
    return compute

def absInt(xInt):
    return Kaucher(abs(xInt.lower),abs(xInt.upper))
def absInt2(xInt):
    return IReal(abs(xInt.inf),abs(xInt.sup))

def diamAbs(xInt):
    return abs(xInt.upper - xInt.lower)

def diamAbs2(xInt):
    return abs(xInt.sup - xInt.inf)


def mKaucher(xInt):
    return (xInt.upper+xInt.lower)/2
def mInt(xInt):
    return (xInt.sup+xInt.inf)/2

def newtons_method(f, x, dx=0.000001, tolerance=0.000000001):
    '''f is the function f(x)'''
    df = derivative(f)
    i = 0.0
    while i<100:
        i+=1.0
        x1 = x - f(x)/df((x), dx)
        t = abs(x1 - x)
        if(t<tolerance):
            break
        x = x1
    return x

def newtons_method_interval(f, x, dx=0.000001, tolerance=0.000000001):
    '''f is the function f(x)'''
    df = derivative(f)
    i = 0.0
    t=0.0
    x1=0.0
    while i<100:
        i+=1.0
        if(type(x) is Kaucher):
            x1 = x - f(x)/df(mKaucher(x), dx)
            t = absInt(x1 - x)
        else:
            x1 = x - f(x)/df(x, dx)
            t = absInt2(x1 - x)
            t = mInt(t)
        if(t<tolerance):
            break
        x = x1
    return x
def f2(x):
    '''
    here solve x for ...
    3*x**5 - 2*x**3 + 1*x - 37 = 0
    '''
    return powI(x,5)*3 - powI(x,3)*2 + x - IReal(37)

def f(x):
    '''
    here solve x for ...
    3*x**5 - 2*x**3 + 1*x - 37 = 0
    '''
    return x**5*3 - x**3*2 + x - 37


x_approx = 1.0
# f refers to the function f(x)
x_real = newtons_method(f, 1.0)
x_kaucher = newtons_method_interval(f, Kaucher(1.0))
x_int = newtons_method_interval(f2, IReal(1.0))

print("Solve for x in 3*x**5 - 2*x**3 + 1*x - 37 = 0")
print x_real
print x_kaucher
print x_int
''' result ...
Solve for x in 3*x**5 - 2*x**3 + 1*x - 37 = 0
x = 1.722575335786
'''
# optional test (result should be close to zero)
# change dx and tolerance level to make it a little closer
print("Testing with the above x value ...")
