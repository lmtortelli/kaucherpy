''' Newton_method3.py
find value for x in f(x) = 0 using Newton's method
see:
http://en.wikipedia.org/wiki/Newton%27s_method
tested with Python27, IronPython273 and Python33 by vegaseat 08jan2013
'''
from kaucherpy import *
from kaucherpy.core import *
from kaucherpy.kaucher import *
from kaucherpy.utils import QualitativeMetrics as qm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np


def plot():
    plt.switch_backend('agg')

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    # Make data.
    X = np.arange(-100, 100, 10.0)
    Y = np.arange(-5, 5, 1.0)
    X, Y = np.meshgrid(X, Y)
    R = 2*X**3 - 2.5*X - 5
    Z = R

    # Plot the surface.

    # Customize the z axis.
    #ax.plot_wireframe(X, Y, Z, rstride=25, cstride=10)
    # Plot the surface
    ax.plot_surface(X, Y, Z, color='b')
    # Add a color bar which maps values to colors.


    plt.savefig('foo5.png')

#from intpy import *

#def powI(xInt,x):
#    return IReal(xInt.inf**x,xInt.sup**x)

def derivative(f):
    def compute(x, dx):
        #print (x)
        #print ((f(x+dx) - f(x))/dx)
        #print (0 in Kaucher(dx)) 
        #print (f(x+dx) - f(x))
        #print (Kaucher(dx))
        r = (f(x+dx) - f(x))/Kaucher(dx)
        print ("OK "+str(r))
        return r
    return compute

def absInt(xInt):
    return Kaucher(abs(xInt.lower),abs(xInt.upper))

def diamAbs(xInt):
    return abs(xInt.upper - xInt.lower)


def mKaucher(xInt):
    return (xInt.upper+xInt.lower)/2

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

def newtons_method_interval(f, x, dx=0.01, tolerance=0.01):
    '''f is the function f(x)'''
    df = derivative(f)
    i = 0.0
    t=0.0
    x1=0.0
    while i<100:
        i+=1.0
        if(type(x) is Kaucher):
            x1 = x - f(x)/df(x, dx)
            t = absInt(x1 - x)
        if(t<tolerance):
            x=x1
            break
        x = x1
    return x

def f(x):
    '''
    here solve x for ...
    3*x**5 - 2*x**3 + 1*x - 37 = 0
    '''
    return 2*x**3 - 2.5*x - 5

plot()
x_approx = 1.0
# f refers to the function f(x)
#x_real = newtons_method(f, 1.0)
x_kaucher = newtons_method_interval(f, Kaucher(8,9.0))
#x_int = newtons_method_interval(f2, IReal(1.0))



print("Solve for x in 3*x**5 - 2*x**3 + 1*x - 37 = 0")
#print x_real
#print (f)
print (x_kaucher)
#print x_int
''' result ...
Solve for x in 3*x**5 - 2*x**3 + 1*x - 37 = 0
x = 1.722575335786
'''
# optional test (result should be close to zero)
# change dx and tolerance level to make it a little closer
print("Testing with the above x value ...")
