import io, subprocess, sys  # importo delle librerie core di python
from math import e 


try:
    
    import numpy as np 
    from sympy import *
    from scipy.optimize import minimize_scalar
    
except:
    
    subprocess.call([sys.executable, "-m", "pip", "install", "numpy", "--break-system-packages"])
    subprocess.call([sys.executable, "-m", "pip", "install", "sympy", "--break-system-packages"])
    subprocess.call([sys.executable, "-m", "pip", "install", "scipy", "--break-system-packages"])
    
    import numpy as np
    from sympy import *
    from scipy.optimize import minimize_scalar
    
    


def Lagrangian(z: float, x: list[int], y: list[float]) -> float :
    
    px = 0
    
    for i in range(len(y)):
        
        c = y[i]
        
        for j in range(len(x)):
            
            if i != j:
            
                c *= (z - x[j]) / (x[i] - x[j])
                
        px += c
    
    return px




def Newtonian(z: float, x: list[int], y: list[float]) -> float :
    
    return


def E(x: float, f, s, ab: list[float]) -> float:
    """
    Calcola una stima dell'errore di interpolazione usando la formula:
    |f(x) - P(x)| ≤ |f^(n+1)(ξ)|/(n+1)! * |w(x)|
    dove w(x) = ∏(x - x_i) e ξ ∈ [min(x_i, x), max(x_i, x)]
    """
    n = len(ab) - 1  # grado del polinomio
    
    # Calcola w(x)
    w_x = 1.0
    for x_i in ab:
        w_x *= (x - x_i)
    
    # Calcola la derivata (n+1)-esima di f
    f_deriv = f
    for i in range(n + 1):
        f_deriv = f_deriv.diff(s)
    
    f_deriv_numeric = lambdify(s, f_deriv, 'numpy')
    
    # Trova il massimo di |f^(n+1)| nell'intervallo appropriato
    interval = [min(min(ab), x), max(max(ab), x)]
    
    # Trova il massimo della derivata nell'intervallo
    result = minimize_scalar(lambda t: -abs_deriv(t, f_deriv, s), bounds=interval, method='bounded')
    max_deriv = abs_deriv(result.x, f_deriv, s)
    
    # Calcola l'errore
    error_bound = (max_deriv / factorial(n + 1)) * abs(w_x)
    
    return error_bound


def abs_deriv(t, f, s) -> float:

    f = lambdify(s, f, "numpy")
    
    return abs(float(f(t)))

def tronca(n: float, i: int) -> float:
    return float(f"{n:.{i}f}")

        
        
def Horner(a: list[float], x: float) -> float:
    
    p = a[-1]
    
    for i in range(len(a)-1):
        
        p = (p * x) + a[((len(a)-2) - i)]
    
    return p


def Taylor(f, s, x: float, x0: float, n: int) -> float:
    
    tmp = lambdify(s, f, 'numpy')
    der = f.diff(s)
    d = lambdify(s, der, 'numpy')
    
    t = tmp(x0)
    
    
    for i in range(1, n+1):
        
        t += d(x0) * (((x-x0)**(i))/factorial((i)))
        
        der = der.diff(s)
        d = lambdify(s, der, 'numpy')
    
    return t


def LinearSpline(x: float, xn: tuple[float, float], yn: tuple[float, float]) -> float:
    
    y = yn[0] + (x-xn[0])*((yn[0]/(xn[0]-xn[1]) + (yn[1]/(xn[1]-xn[0]))))
    
    return y


def FindeMinNSubIntervalls(interval: tuple[float, float], f, s, t: int=(10**-3), max: int=10000, p: int=5) -> tuple[float, float, int]:
    
    h = -1.0
    err = -1.0
    
    for i in range(1, max, p):
        
        dec = np.linspace(interval[0], interval[1], i)
        
        err = -1
        
        for j in range(1, len(dec)):
            
            if E((tronca((dec[j]+dec[j-1])/2, 15)), f, s, [tronca(dec[j-1], 15), tronca(dec[j], 15)]) > err:
                
                err = E((tronca((dec[j]+dec[j-1])/2, 15)), f, s, [tronca(dec[j-1], 15), tronca(dec[j], 15)])
                h = (dec[-1]-dec[0])/(len(dec)-1)
            
                        
        if err <= t and err != -1:
            
            return (tronca(h, int(len(dec)/10)), err, i)
            
    h = -1.0
    err = -1.0
    
    return (h, err, -1)


def binomial(k: int, n: int) -> float:
    
    return (factorial(n) / (factorial(k) * factorial(n-k)))


def bernstein(k: int, n: int, x: float) -> float :
    
    return (binomial(k, n) * (x**k) * ((1-x)**(n-k)))




def L(vs: list[float], v: float, i: int) -> float:
    
    l = 1
    
    for j in range(len(vs)):
        
        if j != i:
            
            l *= ((v-vs[j]) / (vs[i]-vs[j]))
            
    return l




def Lagrangian_3D(xs: list[float], ys: list[float], fxy, s, t, x: float, y: float) -> float :  # funzione che calcola il polinomio di Lagrange in R2
    
    pxy = 0
    
    fxy = lambdify([s, t], fxy, "numpy")
    
    
    for i in range(len(xs)):
        
        d = 0
        
        for j in range(len(ys)):
            
            c = fxy(xs[i], ys[j])
            
            d += c * L(ys, y, j) 
        
        pxy += d * L(xs, x, i) 
                
    
    return pxy