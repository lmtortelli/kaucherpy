#!/usr/bin/env python

    # zillions of papers and methods for derivative-free optimization alone


#...............................................................................
from __future__ import division
import numpy as np
from numpy import abs, cos, exp, mean, pi, prod, sin, sqrt, sum
from intpy import *
import time
import math
#try:
#    from opt.testfuncs.powellsincos import Powellsincos
#except ImportError:
#    Powellsincos = None
#try:
#    from opt.testfuncs.randomquad import randomquad
#    from opt.testfuncs.logsumexp import logsumexp
#except ImportError:
#    randomquad = logsumexp = None


#__version__ = "2015-03-06 mar  denis-bz-py t-online de"  # + randomquad logsumexp

def degToRad(value):
    rad = (IReal(value.inf/180.0,value.sup/180.0))*math.pi
    if (type(rad) is IReal):
        return rad
    else:
        return IReal(rad)
    

def exp(x):
    return IReal(math.e**x.inf,math.e**x.sup)


def sqrt(x):
    if(type(x) is IReal):
        return IReal(x.inf**(1.0/2.0),x.sup**(1.0/2.0))
    else:
        return IReal(x**(1.0/2.0))

def abs(x):
    if(type(x) is IReal):
        return IReal(+abs(x.inf),+abs(x.sup))
    else:
        return np.abs(x)

def sin(value):
    r = 0
    r = (IReal((-cos(value)).inf**2,(-cos(value)).sup**2) + 1)
    r_1 = IReal(r.inf**(1.0/2.0),r.sup**(1.0/2.0))
    return r_1

def cos(rad):
    n = degToRad(rad)
    c = 0
    r = IReal(1)
    while(c < 50):
        c += 1
        t2 = IReal(n.inf**(2 * c),n.sup**(2 * c))
        t3 = IReal(t2.inf/(math.factorial(2 * c)),t2.sup/(math.factorial(2 * c)))
        r = r + IReal(((-1)**c)*t3.inf,((-1)**c)*t3.sup)
    return r

#...............................................................................
def ackley( x, a=20, b=0.2, c=2*math.pi ):
    #x = np.asarray_chkfinite(x)  # ValueError if any NaN or Inf
    #print (x)
    n = len(x)
    s1 = sum( [IReal(v.inf**2,v.sup**2) for v in x])
    s2 = sum( [cos( v * c ) for v in x])
    t1_1 = IReal(s1.inf / n,s1.sup / n)
    t1_2 = IReal(s2.inf / n,s2.sup / n)
    t1_3 = IReal(t1_1.inf**(1.0/2.0),t1_1.sup**(1.0/2.0))
    t2_1 = IReal(b*( t1_3.inf),b*( t1_3.sup))
    return -IReal(a)*exp( -t2_1) - exp( t1_2 ) + IReal(a) + exp(IReal(1))

#...............................................................................
def dixonprice( x ):  # dp.m
    x = np.asarray_chkfinite(x)
    n = len(x)
    j = np.arange( 2, n+1 )
    x2 = x**2 * 2
    t2 = IReal(0.0)
    for i in range(len(x)-1):
        t2+= x2[i+1] - x[i]
    return sum( IReal(t2.inf**2,t2.sup**2)* j ) + (IReal((x[0].inf - 1)**2,(x[0].sup - 1)**2))

#...............................................................................
def griewank( x, fr=4000 ):
    #x = np.asarray_chkfinite(x)
    n = len(x)
    j = np.arange( 1., n+1 )
    s = sum( [IReal(v.inf**2,v.sup**2) for v in x] )
    v = []
    for i in range(len(x)):
        v.append(x[i]/sqrt(j[i]))
    p = prod(([ cos( xi ) for xi in v ]) ) 
    return s/fr - p + 1

#...............................................................................
'''
def levy( x ):
    x = np.asarray_chkfinite(x)
    n = len(x)
    z_1 = (x - 1) + 1
    z = IReal(z_1.inf/4.0,z_1.sup/4.0)
    t1 = sin( pi * z[0] )**2
    t2_1 = [(v - 1)**2 for v in x]
    t2_1 = t2_1[:-1]
    t2_2 = [(1 + 10 * sin( math.pi * v + 1 )**2 ) for v in x]
    t2_2 = t2_2[:-1]
    t2_3 = [t2_1[i]  * t2_2[i] for i in range(len(t2_2))]

    t2 = sum(t2_3)
    t3 = (z[-1] - 1)**2 * (1 + sin( 2 * math.pi * z[-1] )**2 )
    return t1 + t2 + t3      
'''
#...............................................................................
michalewicz_m = .5  # orig 10: ^20 => underflow

def michalewicz( x ):  # mich.m
    #x = np.asarray_chkfinite(x)
    n = len(x)
    j = np.arange( 1., n+1 )
    t1 = [sin(v) for v in x]
    t2 = [sin(IReal(x[i].inf**2/ math.pi,x[i].sup**2/ math.pi ) * j[i]) for i in range(len(x))]
    t3 = [IReal(t2[i].inf**((2 * michalewicz_m)),t2[i].sup**((2 * michalewicz_m)))*t1[i]  for i in range(len(t1))]

    return - sum(t3)

#...............................................................................
'''
def perm( x, b=.5 ):
    x = np.asarray_chkfinite(x)
    n = len(x)
    j = np.arange( 1., n+1 )
    xbyj =  [((IReal(abs(x[i]).inf / j[i],abs(x[i]).sup / j[i]))) for i in range(len(x))]
    return media([ media( (j**k + b) * IReal((xbyj.inf ** k - 1) **2,(xbyj.sup ** k - 1)  **2)) for k in j/n ])
    # original overflows at n=100 --
    # return sum([ sum( (j**k + b) * ((x / j) ** k - 1) ) **2
    #       for k in j ])
'''
#...............................................................................
def powell( x ):
    x = np.asarray_chkfinite(x)
    n = len(x)
    n4 = ((n + 3) // 4) * 4
    if n < n4:
        x = np.append( x, np.zeros( n4 - n ))
    x = x.reshape(( 4, -1 ))  # 4 rows: x[4i-3] [4i-2] [4i-1] [4i]
    f = np.empty_like( x )
    f[0] = x[0] + x[1] * 10
    f[1] = sqrt(5) * (x[2] - x[3])
    f[2] = (x[1] - 2 * x[2]) **2
    f[3] = sqrt(10) * (x[0] - x[3]) **2
    return sum( f**2 )

#...............................................................................
def powersum( x, b=[8,18,44,114] ):  # power.m
    x = np.asarray_chkfinite(x)
    n = len(x)
    s = IReal(0)
    
    for k in range( 1, n+1 ):
        s_1 = IReal(0.0)
        bk = b[ min( k - 1, len(b) - 1 )]  # ?
        for x_1 in x:
            s_1 += (sum( IReal(x_1.inf**k,x_1.inf**k) ) - bk)   # dim 10 huge, 100 overflows
        s+= IReal(s_1.inf**2,s_1.sup**2)
    return s

#...............................................................................
def rastrigin( x ):  # rast.m
    x = np.asarray_chkfinite(x)
    n = len(x)
    return IReal(10*n) + sum([ IReal(v.inf**2,v.sup**2) - cos( v * 2 * math.pi ) * 10 for v in x  ])

#...............................................................................
'''
def rosenbrock( x ):  # rosen.m
    """ http://en.wikipedia.org/wiki/Rosenbrock_function """
        # a sum of squares, so LevMar (scipy.optimize.leastsq) is pretty good
    x = np.asarray_chkfinite(x)
    x0 = x[:-1]
    x1 = x[1:]
    r = IReal(0.0)
    for k in range (len(x0)):
        t1 = sum( (x1[k] - IReal(x0[k].inf**2,x0[k].sup**2)) **2 )
       r+= (sum( IReal((-x0[k].inf + 1) **2,(-x0[k].sup + 1) **2) + IReal(100) * t1 ))
    
'''

#...............................................................................
def schwefel( x ):  # schw.m
    x = np.asarray_chkfinite(x)
    n = len(x)
    return IReal(418.9829*n) - sum( [v * sin( sqrt( abs( v ))) for v in x] )

#...............................................................................
def sphere( x ):
    x = np.asarray_chkfinite(x)
    return sum( x**2 )

#...............................................................................
def sum2( x ):
    x = np.asarray_chkfinite(x)
    n = len(x)
    j = np.arange( 1., n+1 )
    x_1 = [IReal(v.inf**2,v.sup**2) for v in x]
    return sum( x_1 * j )

#...............................................................................
def trid( x ):
    x = np.asarray_chkfinite(x)
    return sum( (x - 1) **2 ) - sum( x[:-1] * x[1:] )

#...............................................................................
def zakharov( x ):  # zakh.m
    x = np.asarray_chkfinite(x)
    n = len(x)
    j = np.arange( 1., n+1 )
    s2_1 = sum( [IReal(v.inf * j,v.sup * j ) for v in x])
    s2 = IReal(s2_1.inf/ 2,s2_1.sup/ 2)
    return sum( x**2 ) + s2**2 + s2**4

#...............................................................................
    # not in Hedar --
def media(x):
    t = sum(x)
    return IReal(t.inf/len(x),t.sup/len(x))
    
def ellipse( x ):
    x = np.asarray_chkfinite(x)
    t2 = np.diff(x)
    return media( (-x + 1) **2 )  + media( t2 **2 ) * 100
#...............................................................................

def m(x):
    return (x.inf+x.sup)/2.0

def nesterov( x ):
    """ Nesterov's nonsmooth Chebyshev-Rosenbrock function, Overton 2011 variant 2 """
    x = np.asarray_chkfinite(x)
    x0 = x[:-1]
    x1 = x[1:]
    t1 = abs( 1 - m(x[0]) ) / 4
    t2 = sum( [abs(m( x1[i] - 2*abs(m(x0[i])) + 1 )) for i in range(len(x1)) ])
    return t1 + t2

#...............................................................................
def saddle( x ):
    x = np.asarray_chkfinite(x) - 1
    return np.mean( np.diff( x **2 )) \
        + .5 * np.mean( x **4 )


#-------------------------------------------------------------------------------
allfuncs = [
    ackley,
    dixonprice,
    ellipse,
 #   griewank,
  #  levy,
    michalewicz,  # min < 0
    nesterov,
  #  perm,
  #  powell,
    # powellsincos,  # many local mins
    powersum,
    rastrigin,
#    rosenbrock,
    schwefel,  # many local mins
    sphere,
    #saddle,
    sum2,
    trid,  # min < 0
 #   zakharov,
    ]

'''
if Powellsincos is not None:  # try import
    _powellsincos = {}  # dim -> func
    def powellsincos( x ):
        x = np.asarray_chkfinite(x)
        n = len(x)
        if n not in _powellsincos:
            _powellsincos[n] = Powellsincos( dim=n )
        return _powellsincos[n]( x )

    allfuncs.append( powellsincos )
'''

#if randomquad is not None:  # try import
#    allfuncs.append( randomquad )
#    allfuncs.append( logsumexp )

allfuncs.sort( key = lambda f: f.__name__ )


#...............................................................................
allfuncnames = " ".join([ f.__name__ for f in allfuncs ])
name_to_func = { f.__name__ : f  for f in allfuncs }

    # bounds from Hedar, used for starting random_in_box too --
    # getbounds evals ["-dim", "dim"]
ackley._bounds       = [-15, 30]
dixonprice._bounds   = [-10, 10]
griewank._bounds     = [-600, 600]
#levy._bounds         = [-10, 10]
michalewicz._bounds  = [0, pi]
#perm._bounds         = ["-dim", "dim"]  # min at [1 2 .. n]
#powell._bounds       = [-4, 5]  # min at tile [3 -1 0 1]
powersum._bounds     = [0, "dim"]  # 4d min at [1 2 3 4]
rastrigin._bounds    = [-5.12, 5.12]
#rosenbrock._bounds   = [-2.4, 2.4]  # wikipedia
schwefel._bounds     = [-500.0, 500.0]
sphere._bounds       = [-5.12, 5.12]
sum2._bounds         = [-10, 10]
trid._bounds         = ["-dim**2", "dim**2"]  # fmin -50 6d, -200 10d
#zakharov._bounds     = [-5, 10]

ellipse._bounds      =  [-2, 2]
#logsumexp._bounds    = [-20, 20]  # ?
nesterov._bounds     = [-2, 2]
#powellsincos._bounds = [ "-20*pi*dim", "20*pi*dim"]
#randomquad._bounds   = [-10000, 10000]
#saddle._bounds       = [-3, 3]


#...............................................................................
def getfuncs( names, dim=0 ):
    """ for f in getfuncs( "a b ..." ):
            y = f( x )
    """
    if names == "*":
        return allfuncs
    funclist = []
    for nm in names.split():
        if nm not in name_to_func:
            raise ValueError( "getfuncs( \"%s\" ) not found" % names )
        funclist.append( name_to_func[nm] )
    return funclist

def getbounds( funcname, dim ):
    """ "ackley" or ackley -> [-15, 30] """
    funcname = getattr( funcname, "__name__", funcname )
    func = getfuncs( funcname )[0]
    try:
        b = func._bounds[:]
        #print (b)
        if isinstance( b[0], str ):  b[0] = eval( b[0] )
        if isinstance( b[1], str ):  b[1] = eval( b[1] )
        return b
    except:
        return None

    

#...............................................................................
_minus = "dixonprice perm powersum schwefel sphere sum2 trid zakharov "  # nlopt ~ same ?

def allfuncs_minus( minus=_minus ):
    return [f for f in allfuncs
            if f.__name__ not in minus.split() ]

def funcnames_minus( minus=_minus ):
    return " ".join([ f.__name__
            for f in allfuncs_minus( minus=minus )])

def makeList(steps,dim):
    result = []
    r = []
    for k in steps:
        for i in range(dim):
            r.append(IReal(k))
        result.append(r)
        r = []

    return result
#-------------------------------------------------------------------------------
if __name__ == "__main__":  # standalone test --
    import sys
    import json

    dims = [2,4,6,8]
    nstep = 700  # 11: 0 .1 .2 .. 1
    seed = 1
    problems = {}
    init = 50
    
        # to change these params in sh or ipython, run this.py  a=1  b=None  c=[3] ...
    for arg in sys.argv[1:]:
        exec( arg )

    np.set_printoptions( threshold=20, edgeitems=5, linewidth=120, suppress=True,
        formatter = dict( float = lambda x: "%.2g" % x ))  # float arrays %.2g
    np.random.seed(seed)

    #...........................................................................
    for dim in dims:


        print ("\n# ndtestfuncs dim %d  along the diagonal low .. high corner --" % dim)
        # cmp matlab, anyone ?
        print ("PROBLEM | DIM | MIN | LOW | HI | Y - YMIN")
        for step in range(init,nstep):
            print ("-Dims: #",dim,"  -STEP: #",step)
            for func in allfuncs:
                lo, hi = getbounds( func, dim )
                stepsAux = np.linspace( lo, hi, step )
                steps = makeList(stepsAux,dim)
                startTime = time.time()  
                Y = np.array([ func(t) for t in steps ])
                endTime = time.time()
                jmin = Y.argmin()
                Ymin = Y[jmin]
                if(step==init and dim == dims[0]):
                    problems[func.__name__] = {"dim":dim,"Ymin":Ymin,"cpu":(endTime - startTime),"Low":lo,"High":hi,"N-Steps":step}
                else:
                    if(problems[func.__name__]["Ymin"] > Ymin):
                        problems[func.__name__] = {"dim":dim,"Ymin":Ymin,"cpu":(endTime - startTime),"Low":lo,"High":hi,"N-Steps":step}
                #print ("problem dim min lo hi Y-min ",
                #        func.__name__, dim, Ymin, steps[jmin], lo, hi, Y - Ymin )
                #print (func.__name__,dim,Ymin,lo,hi,Y-Ymin)
                #print (problems)
    with open("intervalIntpy.txt", 'w') as f:
        for key, value in problems.items():
            f.write('%s:%s\n' % (key, value))
    f.close()