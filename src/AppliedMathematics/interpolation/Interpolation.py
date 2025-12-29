import numpy as np 
from sympy import *
from sympy.core.function import UndefinedFunction
from scipy import *
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core.function import UndefinedFunction
from numpy import float64

from math import factorial


def LS(x: float64, xn: list[float64], yn: list[float64]) -> float64:
    
    y = float64(yn[0] + (x-xn[0])*((yn[0]/(xn[0]-xn[1]) + (yn[1]/(xn[1]-xn[0])))))
    
    return y


def Linear_Spline(x: float64, xn: list[float64], yn: list[float64]) -> float64:
    
    if len(xn) != len(yn):
        
        raise Exception("Sorry, xn and yn must have the same lenght")
    
    if x > max(xn):
        
        raise Exception("Sorry x is not in the range")
    
    idx = np.searchsorted(xn, x) - 1
    
    return LS(x, xn[idx:idx+2], yn[idx:idx+2])


def Symbolic_Linear_Spline(s: Symbol, xn: list[float64], yn: list[float64]) -> UndefinedFunction:
    
    if len(xn) != len(yn):
        
        raise Exception("Sorry, xn and yn must have the same lenght")
    
    
    n = len(xn)
    
    sp = 0
    
    for i in range(n-1):
        
        sp = yn[i] + (s-xn[i])*((yn[i]/(xn[i]-xn[i+1]) + (yn[i+1]/(xn[i+1]-xn[i]))))
        
    
    return sp




def Taylor(f, s, x: float64, x0: float64, points: int) -> float64:
    
    n = points
    
    tmp = lambdify(s, f, 'numpy')
    der = f.diff(s)
    d = lambdify(s, der, 'numpy')
    
    t = tmp(x0)
    
    
    for i in range(1, n+1):
        
        t += float64(d(x0) * (((x-x0)**(i))/factorial((i))))
        
        der = der.diff(s)
        d = lambdify(s, der, 'numpy')
    
    return t


def Symbolic_Taylor(f, s, x: Symbol, x0: float64, points: int) -> Add:
    
    
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



def L(vs: list[float64], v: float64, i: int) -> float64:
    
    l = float64(1)
    
    for j in range(len(vs)):
        
        if j != i:
            
            l *= float64((v-vs[j]) / (vs[i]-vs[j]))
            
    return l


def Lagrange(x: float64, xn: list[float64], yn: list[float64]) -> float64:
    
    if len(xn) != len(yn):
        
        raise Exception("Sorry, xn and yn must have the same lenght")
    
    
    l = float64(0)
    
    n = len(xn)
    
    for i in range(n):
        
        l += float64(yn[i] * L(xn, x, i))
        

    return l


def Symbolic_L(vs: list[float64], s: Symbol, i: int) -> Add:
    
    l = float64(1)
    
    for j in range(len(vs)):
        
        if j != i:
            
            l *= ((s-vs[j]) / (vs[i]-vs[j]))
            
    return l



def Symbolic_Lagrange(s: Symbol, xn: list[float64], yn: list[float64]) -> Add:
    
    
    if len(xn) != len(yn):
        
        raise Exception("Sorry, xn and yn must have the same lenght")
    
    
    l = float64(0)
    
    n = len(xn)
    
    for i in range(n):
        
        l += yn[i] * Symbolic_L(xn, s, i)
        

    return l



def Divided_Differences(x: list[float64], y: list[float64]) -> float64:
    
    f = float64(0)
    
    n = len(x)
    
    for i in range(n):
        
        m = 1
        
        for j in range(n):
            
            if i != j:
                
                m *= float64(x[i] - x[j])
                
        f += float64(y[i] / m)
    
    return f



def Symbolic_Newtonian_Polynomials(x: list[float64], y: list[float64], s: Symbol) -> Add:
    
    if len(x) != len(y):
        
        raise Exception("Sorry, x and y must have the same lenght")
    
    n = len(x)
    
    p = Divided_Differences([x[0]], [y[0]])
    
    for i in range(1, n):
        
        tmp = 1
        
        for j in range(i):
            
            tmp *= (s-x[j])
            
        p += tmp*Divided_Differences(x[:i+1], y[:i+1])
        
    return p


def Newtonian_Polynomials(z: float64, x: list[float64], y: list[float64]) -> float64:
    
    if len(x) != len(y):
        
        raise Exception("Sorry, x and y must have the same lenght")
    
    n = len(x)
    
    p = Divided_Differences([x[0]], [y[0]])
    
    for i in range(1, n):
        
        tmp = 1
        
        for j in range(i):
            
            tmp *= float64(z-x[j])
            
        p += float64(tmp*Divided_Differences(x[:i+1], y[:i+1]))
        
    return p



# ToDo errore spline, spline cubiche, quadratiche
# ToDo Gershgorin
# ToDo interpolazione a due varibili con spline



if __name__ == "__main__":
    
    y = Symbol("y")
    z = Symbol("z")
    
    f = y
    
    x = float64(3.0)
    
    xn = [0.0, 1, 2, 4, 5]
    yn = [0.0, 1, 2, 4, 5]
    
    
    
    print("#"*10, end="")
    print(" Info ", end="")
    print("#"*10)
    
    print(f"Symbols: {y}, {z}")
    print(f"Function: f(y)={f}")
    print(f"Number to find: {x}")
    print(f"Points: {[(xn[i], yn[i]) for i in range(len(xn))]}", end="\n\n\n")
    
    
    print("#"*10, end="")
    print(" Linear Spline and Symbolic Linear Spline", end="")
    print("#"*10)
    
    print(f"Lineal Spline: {Linear_Spline(x, xn, yn)}")
    sp = Symbolic_Linear_Spline(y, xn, yn)
    sp = lambdify(y, sp, "numpy")
    print(f"Symbolic Lineal Spline: {sp(x)}", end="\n\n\n")
    
    
    print("#"*10, end="")
    print(" Taylor and Symbolic Taylor ", end="")
    print("#"*10)
    
    print(f"Taylor: {Taylor(f, y, x, xn[0], len(xn))}")
    
    t = Symbolic_Taylor(f, y, z, xn[0], len(xn))
    
    t = lambdify(z, t, "numpy")
    
    print(f"Symbolic Taylor: {t(x)}", end="\n\n\n")
    
    
    print("#"*10, end="")
    print(" Lagrange and Symbolic Lagrange ", end="")
    print("#"*10)
    
    
    print(f"Lagrange: {Lagrange(x, xn, yn)}")
    
    l = Symbolic_Lagrange(y, xn, yn)
    
    l = lambdify(y, l, "numpy")
    
    print(f"Symbolic Lagrange: {l(x)}", end="\n\n\n")
    
    
    print("#"*10, end="")
    print(" Newtonian and Symbolic Newtonian ", end="")
    print("#"*10)
    
    np = Symbolic_Newtonian_Polynomials(xn, yn, y)
    
    
    np = lambdify(y, np, "numpy")
    
    print(f"Newtonian: {Newtonian_Polynomials(x, xn, yn)}")
    
    print(f"Symbolic Newtonian: {np(x)}")
    
    
    