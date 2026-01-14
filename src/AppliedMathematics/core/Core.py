from math import factorial
from numpy import float64
from sympy import Symbol
from sympy.core.function import UndefinedFunction


nullsym = Symbol("")
nullfunc = nullsym
nullfunc = nullfunc.diff(nullsym, 2)
nullfloat = float64(None)


def Horner(a: list[float64], x: float64) -> float64:
    
    p = a[-1]
    
    for i in range(len(a)-1):
        
        p = (p * x) + a[((len(a)-2) - i)]
    
    return p



def binomial(k: int, n: int) -> float64:
    
    return float64(factorial(n) / (factorial(k) * factorial(n-k)))


def convert(x: float64) -> tuple[int, float64]:
    
    i = 0
    
    while x < 1:
        
        x = float64(x * 10)
        
        i += 1
    
    return (i, x)


def tronca(n: float64, i: int) -> float64:
    
    return float64(f"{n:.{i}f}")



if __name__ == "__main__":
    
    a = [1, 2, 3, 7.0]
    
    x = float64(2)
    
    print(Horner(a, x))