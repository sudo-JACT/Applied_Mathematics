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
    
    
    
def wi(n: int) -> list[float]:  # funzione che restituisce i pesi di Newton-Cotes
    
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



def Newton_Cotes(n: int, f, s: Symbol, a: float, b: float) -> float:  # fuzione che calcola la quadratura di Newton-Cotes
    
    f = lambdify(s, f, "numpy")
    
    y = f(np.linspace(a, b, n+1))  # Calcolo tutti gli yi
    
    w = wi(n)
    
    h = (b-a) / n
    
    I = 0
    
    for i in range(n+1):
        
        I += (y[i] * w[i])
        
    
    return I * h