import numpy as np 
from sympy import *
from sympy.core.function import UndefinedFunction
from scipy import *




def Cavalieri_Simpson(a: float, b: float, f, s: Symbol) -> float:
    
    f = lambdify(s, f, "numpy")
    
    m = (a+b) / 2
    
    return (((b-a)/6) * (f(a) + 4*f(m) + f(b)))
    



def Cavalieri_Simpson_Comp(a: float, b: float, ep: float, n: int, f, s: Symbol) -> tuple[int, float, float]:

    
    fIV = f.diff(s, 4)
    
    f1 = lambdify(s, f, "numpy")
    fIV = lambdify(s, fIV, "numpy")

    
    max = np.max(np.abs(fIV(np.linspace(a, b, 10000))))
    
    
    E = abs((((b-a)**5) / (2880 * ((n)**4))) * max)
    
    
    if E <= ep:
        
        x = np.linspace(a, b, n+1)
        
        sum1, sum2 = 0, 0
        
        for i in range(1, n):
            
            sum1 += f1(x[i])
            
        
        for i in range(n):
            
            sum2 += f1((x[i] + x[i+1]) / 2)
        
        
        abn = (b-a) / (6*(n))
        
        fs = f1(a) + (2*sum1) + (4*sum2) + f1(b)
        
        Q = (abn * (fs))
        
    else:
        
        n, Q, E = Cavalieri_Simpson_Comp(a, b, ep, (1+n), f, s)
    
    return (n, Q, E)
    
    
    
def Newton_Cotes_Weights(n: int) -> list[float]:  # funzione che restituisce i pesi di Newton-Cotes
    
    w = [
        [1/2, 1/2],
        [1/3, 4/3, 1/3],
        [3/8, 9/8, 9/8, 3/8],
        [14/45, 64/45, 24/45, 64/45, 14/45],
        [95/288, 375/288, 250/288, 250/288, 375/288, 95/288],
        [41/140, 216/140, 27/140, 272/140, 27/140, 216/140, 41/140],
        [5257/17280, 25039/17280, 9261/17280, 20923/17280, 20923/17280, 9261/17280, 25039/17280, 5257/17280]
    ]
    
    return w[n-1]



def Newton_Cotes(points: int, f, s: Symbol, a: float, b: float) -> float:  # fuzione che calcola la quadratura di Newton-Cotes
    
    n = points
    
    f = lambdify(s, f, "numpy")
    
    y = f(np.linspace(a, b, n+1))  # Calcolo tutti gli yi
    
    w = Newton_Cotes_Weights(n)
    
    h = (b-a) / n
    
    I = 0
    
    for i in range(n+1):
        
        I += (y[i] * w[i])
        
    
    return I * h


if __name__ == "__main__":
    
    x = Symbol("x")
    
    f = x**2
    
    a, b = 0, 1
    
    e = 10**(-3)
    
    n = 3
    
    print(f"Newton-Cotes of f(x)={f}: {Newton_Cotes(n, f, x, a, b)}")
    
    print(f"Cavalieri-Simpson of f(x)={f}: {Cavalieri_Simpson(a, b, f, x)}")
    
    print(f"Iterated Cavalieri-Simpson of f(x)={f}: {Cavalieri_Simpson_Comp(a, b, e, n, f, x)[1]}")