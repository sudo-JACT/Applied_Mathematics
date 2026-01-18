from sympy import Symbol, lambdify
from sympy.core.function import UndefinedFunction
from numpy import float64
import numpy as np


from ..core import nullfloat


def Bisection_Method(a: float64, b:float64, f: UndefinedFunction, s: Symbol, err: float64=float64(0), iter: int=0) -> tuple[float64, int]:
    
    g = lambdify(s, f, "numpy")
    
    
    if (g(a) - float64(0) <= 0) and (g(b) - float64(0) >= 0):
    
        m = (b+a)/float64(2)
        
        if abs(g(m)) <= err:
            
            return (m, iter)
        
        elif (g(m)) < 0:
            
            return Bisection_Method(m, b, f, s, err, iter+1)
        
        elif (g(m)) > 0:
            
            return Bisection_Method(a, m, f, s, err, iter+1)
        
    
        
    return (nullfloat, iter)





def Secant_Method(a: float64, b:float64, f: UndefinedFunction, s: Symbol, err: float64=float64(0)) -> tuple[float64, int]:
    
    g = lambdify(s, f, "numpy")
    
    
    if (g(a) * g(b)) == float64(0):
        
        raise Exception("One of the roots is one of them")
    
    
    xsi = []
    yi = []
    
    xi = g(b)
    
    xsi.insert(0, a)
    yi.insert(0, g(a))
    
    xsi.insert(0, b)
    yi.insert(0, xi)
    
    
    iter = 0
    
    while abs(yi[0]) > err:
        
        m = (yi[0] - yi[1]) / (xsi[0] - xsi[1])        
        
        xi = xsi[0] - (yi[0] / m)
        
        xsi.pop()
        yi.pop()
        
        xsi.insert(0, xi)
        yi.insert(0, g(xi))
        
        iter += 1
        
        
        
    return (xsi[0], iter)
    
    
    
    
def Chord_Method(a: float64, b:float64, f: UndefinedFunction, s: Symbol, err: float64=float64(0), maxiter: int=1000) -> tuple[float64, int]:
    
    g = lambdify(s, f, "numpy")
    
    if (g(a) * g(b)) == float64(0):
        
        raise Exception("One of the roots is one of them")
    
    m = (g(b) - g(a)) / (b - a)
    
    xi = (b+a) / 2
    
    it = 0
    
    
    while abs(g(xi)) > err and it < maxiter:
        
        
        if abs((g(xi) / m)) < err:
            
            break
        
        
        xi = xi - (g(xi) / m)
        
        it += 1
        
        
    return (xi, it)

    
    
    
    
def Newtons_Method(a: float64, b:float64, f: UndefinedFunction, s: Symbol, x0: float64=nullfloat, err: float64=float64(0), maxiter: int=1000) -> tuple[float64, int]:    
    
    
    g = lambdify(s, f, "numpy")
    
    gI = lambdify(s, f.diff(s), "numpy")
    
    
    if (g(a) * g(b)) == float64(0):
        
        raise Exception("One of the roots is one of them")
    

    
    
    if np.isnan(x0):
        
        x0 = (a+b) / float64(2)
    


    xi = x0
    
    
    it = 0

    
    while abs(g(xi)) > err and it < maxiter:

        
        m = gI(xi)
         
        if (abs(m) < (1e-15)):
            
            if abs(g(xi)) < (1e-06):
                
                if a <= xi <= b:
                    
                    return xi
                
            else:
                    
                xi = xi - (np.sign(g(xi)) * (1e-06))
                        
        else:
        
            xi = xi - (g(xi) / m)
    
        

        it += 1
    
    
    
    return (xi, it)




def Regula_Falsi(a: float64, b:float64, f: UndefinedFunction, s: Symbol, x0: float64=nullfloat, err: float64=float64(0), maxiter: int=1000) -> tuple[float64, int]:
    
    
    g = lambdify(s, f, "numpy")
    
    
    if (g(a) * g(b)) == float64(0):
        
        raise Exception("One of the roots is one of them")
    
    
    xi = (b+a) / 2
    
    m = (a*g(xi) - xi*g(a)) / (xi - a)
    
    
    it = 0
    
    
    while abs(g(xi)) > err and it < maxiter:
        
        xi = xi - (g(xi) / m)
        
        m = (a*g(xi) - xi*g(a)) / (xi - a)
        
        it += 1
        
    return (xi, it)



if __name__ == "__main__":
    
    x = Symbol("x")
    
    f = (x**3) - 5
    
    a, b = float64(1), float64(2)
    
    err = float64(1e-10)
    
    
    print(Bisection_Method(a, b, f, x, err, 0), end="\n\n\n")
    
    
    print(Secant_Method(a, b, f, x, err), end="\n\n\n")
    
    
    print(Chord_Method(a, b, f, x, err), end="\n\n\n")
    
    
    print(Newtons_Method(a=a, b=b, f=f, s=x, err=err), end="\n\n\n")
    
    
    print(Regula_Falsi(a=a, b=b, f=f, s=x, err=err), end="\n\n\n")