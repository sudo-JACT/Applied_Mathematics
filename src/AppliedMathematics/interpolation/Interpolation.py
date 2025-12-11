import numpy as np 
from sympy import *
from sympy.core.function import UndefinedFunction
from scipy import *

from math import factorial


def LS(x: float, xn: list[float], yn: list[float]) -> float:
    
    y = yn[0] + (x-xn[0])*((yn[0]/(xn[0]-xn[1]) + (yn[1]/(xn[1]-xn[0]))))
    
    return y


def Linear_Spline(x: float, xn: list[float], yn: list[float]) -> float:
    
    if len(xn) != len(yn):
        
        raise Exception("Sorry, xn and yn must have the same lenght")
    
    idx = np.searchsorted(xn, x) - 1
    
    return LS(x, xn[idx:idx+2], yn[idx:idx+2])




def Taylor(f, s, x: float, x0: float, points: int) -> float:
    
    n = points
    
    tmp = lambdify(s, f, 'numpy')
    der = f.diff(s)
    d = lambdify(s, der, 'numpy')
    
    t = tmp(x0)
    
    
    for i in range(1, n+1):
        
        t += d(x0) * (((x-x0)**(i))/factorial((i)))
        
        der = der.diff(s)
        d = lambdify(s, der, 'numpy')
    
    return t



def L(vs: list[float], v: float, i: int) -> float:
    
    l = 1
    
    for j in range(len(vs)):
        
        if j != i:
            
            l *= ((v-vs[j]) / (vs[i]-vs[j]))
            
    return l


def Lagrange(x: float, xn: list[float], yn: list[float]) -> float:
    
    if len(xn) != len(yn):
        
        raise Exception("Sorry, xn and yn must have the same lenght")
    
    
    l = 0
    
    n = len(xn)
    
    for i in range(n):
        
        l += yn[i] * L(xn, x, i)
        

    return l



# ToDo interpolazione di Newton e Differenze divise
# ToDo errore spline, spline cubiche, quadratiche
# ToDo Gershgorin
# ToDo interpolazione a due varibili con spline



if __name__ == "__main__":
    
    y = Symbol("y")
    
    f = y
    
    x = 3
    
    xn = [0.0, 1, 2, 4, 5]
    yn = [0.0, 1, 2, 4, 5]
    
    
    print(Linear_Spline(x, xn, yn))
    
    print(Taylor(f, y, 3, 0, 5))
    
    print(Lagrange(x, xn, yn))
    
    
    
