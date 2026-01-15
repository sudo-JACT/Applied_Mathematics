from sympy import Symbol, lambdify
from sympy.core.function import UndefinedFunction
from numpy import float64


from ..core import nullfloat


def Bisection_Method(a: float64, b:float64, f: UndefinedFunction, s: Symbol, x: float64=float64(0), err: float64=float64(0), iter: int=0) -> tuple[float64, int]:
    
    g = lambdify(s, f, "numpy")
    
    
    if (g(a) - x <= 0) and (g(b) - x >= 0):
    
        m = (b+a)/float64(2)
        
        if abs(g(m) - x) <= err:
            
            return (m, iter)
        
        elif (g(m) - x) < 0:
            
            return Bisection_Method(m, b, f, s, x, err, iter+1)
        
        elif (g(m) - x) > 0:
            
            return Bisection_Method(a, m, f, s, x, err, iter+1)
        
    
        
    return (nullfloat, iter)