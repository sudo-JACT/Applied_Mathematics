import numpy as np 
from sympy import *
from sympy.core.function import UndefinedFunction
from numpy import float64


def tronca(n: float64, i: int) -> float64:
    
    return float64(f"{n:.{i}f}")


def convert(x: float64) -> tuple[int, float64]:
    
    i = 0
    
    while x < 1:
        
        x = x * 10
        
        i += 1
    
    return (i, x)


def Cavalieri_Simpson(a: float64, b: float64, f, s: Symbol) -> float64:
    
    f = lambdify(s, f, "numpy")
    
    m = (a+b) / 2
    
    return (((b-a)/6) * (f(a) + 4*f(m) + f(b)))
    



def Iterated_Cavalieri_Simpson(a: float64, b: float64, err: float64, points: int, f, s: Symbol) -> tuple[int, float64, float64]:

    n = points
    
    fIV = f.diff(s, 4)
    
    f1 = lambdify(s, f, "numpy")
    fIV = lambdify(s, fIV, "numpy")

    
    max = np.max(np.abs(fIV(np.linspace(a, b, 10000))))
    
    
    E = abs((((b-a)**5) / (2880 * ((n+1)**4))) * max)
    
    
    if E <= err:
        
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
        
        n, Q, E = Iterated_Cavalieri_Simpson(a, b, err, (1+n), f, s)
    
    return (n, Q, E)
    
    

def Iterated_Cavalieri_Simpson_With_Points(points: list[float64], h: float64) -> float64:

    n = len(points)
    
    sum1, sum2 = 0, 0
    
    for i in range(1, n-1):
        
        if i % 2 == 0:
            
            sum2 += points[i]
            
        else:
        
            sum1 += points[i]
    
    
    
    abn = h / 3
    
    fs = points[0] + (2*sum1) + (4*sum2) + points[-1]
    
    Q = (abn * (fs))
        
    
    return Q
    

    
    
def Newton_Cotes_Weights(n: int) -> list[float64]:  # funzione che restituisce i pesi di Newton-Cotes
    
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



def Newton_Cotes(points: int, f, s: Symbol, a: float64, b: float64) -> float64:  # fuzione che calcola la quadratura di Newton-Cotes
    
    n = points
    
    f = lambdify(s, f, "numpy")
    
    y = np.linspace(a, b, n+1)  # Calcolo tutti gli yi
    
    y2 = [f(i) for i in y]
    
    y = y2
    
    w = Newton_Cotes_Weights(n)
    
    h = (b-a) / n
    
    I = 0
    
    for i in range(n+1):
        
        I += (y[i] * w[i])
        
    
    return I * h



def Trapezoidal_Rule(a: float64, b: float64, f, s: Symbol) -> float64:
    
    g = lambdify(s, f, "numpy")
    
    
    return ((b-a) / 2) * (g(a) + g(b))




def Itereted_Trapezoidal_Rule(a: float64, b: float64, points: int, err: float64, f, s: Symbol) -> tuple[float64, float64, int]:
    
    E = Iterated_Trapezoidal_Rule_Error(a, b, points, f, s)
    
    if E <= err:
        
        f = lambdify(s, f, "numpy")
    
        x = [f(i) for i in np.linspace(a, b, points+1)]
        
        h = ((b-a) / (points)) / 2
            
        Tn = (x[0] + x[-1]) + 2*np.sum(x[1:-1])   
            
        return (Tn * h, E, points)
    
    else:
        
        return Itereted_Trapezoidal_Rule(a, b, 2*points, err, f, s)
    



def Cavalieri_Simpson_Error(a: float64, b: float64, f, s: Symbol) -> float64:
    
    fIV = f.diff(s, 4)
    fIV = lambdify(s, fIV, "numpy")

    max = np.max(np.abs(fIV(np.linspace(a, b, 10000))))
    
    E = abs((((b-a)**5) / 180) * max)
    
    return E



def Iterated_Cavalieri_Simpson_Error(a: float64, b: float64, points: int, f, s: Symbol) -> float64:
    
    n = points
    
    fIV = f.diff(s, 4)
    fIV = lambdify(s, fIV, "numpy")

    max = np.max(np.abs(fIV(np.linspace(a, b, 10000))))
    
    E = abs((((b-a)**5) / (2880 * ((n+1)**4))) * max)
    
    return E




def Trapezoidal_Rule_Error(a: float64, b: float64, f, s: Symbol) -> float64:
    
    fII = f.diff(s, 2)
    fII = lambdify(s, fII, "numpy")
    
    h = (b-a)
    
    max = np.max(np.abs(fII(np.linspace(a, b, 10000))))
    
    E = abs(((h**3) / 12) * max)
    
    return E



def Iterated_Trapezoidal_Rule_Error(a: float64, b: float64, points: int, f, s: Symbol) -> float64:
    
    fII = f.diff(s, 2)
    fII = lambdify(s, fII, "numpy")
    
    h = (b-a) / (points+1)
    
    max = np.max(np.abs(fII(np.linspace(a, b, 10000))))
    
    E = abs(((((b-a)**3) / (12 * ((points+1)**2)))) * max)
    
    return E




def Filon_Error(a: float64, b: float64, h: float64, w: float64, f, s: Symbol) -> float64:
    
    
    g = f.diff(s, 4)
    
    g = lambdify(s, g, "numpy")
    
    D = np.linspace(a, b, 10000)
    
    ee = [float64(g(i)) for i in D]
    
    ee = max(ee)
    
    th = float64(w*h)
    
    
    if ee == 0:
        
        r = float64((h**3) * ((b - a) / 12) * (1 - (1 / (16*np.cos(th/4)))) * np.sin(th/2))
    
    else:
    
        r = float64((h**3) * ((b - a) / 12) * (1 - (1 / (16*np.cos(th/4)))) * np.sin(th/2) * ee)
    
    return abs(r)




def Filon(a: float64, b: float64, points: int, error: float64, tol: float64, w: float64, f, s: Symbol, cos: bool=True) -> tuple[float64, float64, float64, int]:
    
    
    h = float64((b-a) / (2*points))
    
    E = Filon_Error(a, b, h, w, f, s)
    
    theta = float64(tronca(w*h, int(convert(tol)[0])))
    
    if tol < float64(10**(-6)):
        
        tol = float64(10**(-6))


    
    if E <= error and theta <= tol:
        
        D = np.arange(a, b+h, h)
        
        #D = [a+i*h for i in range(points+1)]
    
        
        if abs(theta) < 1e-10:
            
            alpha = float64(0)
            beta = float64(1/3)
            gamma = float64(4/3)
            
        else:
            
            alpha = float64(1/theta) + (np.sin(2*theta)/(2*(theta**2))) - ((2*(np.sin(theta))**2))/(theta**3)
            beta = float64((1+(np.cos(theta)**2))/(theta**2)) - (np.sin(2*theta)/(theta**3))
            gamma = float64(4/(theta**3)) * (np.sin(theta) - theta*np.cos(theta))
        
        
        
        
        f = lambdify(s, f, "numpy")
        
        
        if cos:
            
            c2 = 0
            
            for i in range(1, points):
                
                c2 += float64(f(a+(2*i*h))*np.cos(w*(a+(2*i*h))))
            
            C2 = float64((f(a)*np.cos(w*a)) + 2*c2 + f(b)*np.cos(w*b))
            
            C1 = 0
            
            for i in range(1, points+1):
                
                C1 += float64(f(a+h*(2*i-1))*np.cos(w*a + w*h*(2*i-1)))
                
            
            Q = float64((h*alpha * (f(b)*np.sin(w*b) - f(a)*np.sin(w*a))) + beta*h*C2 + gamma*h*C1)
                
        else:
            
            s2 = 0
            
            for i in range(1, points):
                
                s2 += float64(f(a+2*i*h)*np.sin(w*(a+2*i*h)))
            
            S2 = float64((f(a)*np.sin(w*a)) + s2 + f(b)*np.sin(w*b))
            
            S1 = 0
            
            for i in range(1, points+1):
                
                S1 += float64(f(a+h*(2*i-1))*np.sin(w*h*(2*i-1)))
                
            
            Q = float64((h*alpha * (f(a)*np.cos(w*a) - f(b)*np.cos(w*b))) + beta*h*S2 + gamma*h*S1)
        
        
        
        return (Q, E, h, points)
        
        
    
    else:
        
        return Filon(a, b, 10*points, error, tol, w, f, s, cos)
    
    
    
    






if __name__ == "__main__":
    
    x = Symbol("x")
    
    f = x**2
    
    #f = f.diff(x)
    
    a, b = float64(0), float64(1)
    
    e = float64(10**(-8))
    
    n = 3
    
    
    print(f"Cavalieri-Simpson Error of f(x)={f}: {Cavalieri_Simpson_Error(a, b, f, x)}", end="\n\n")
    
    print(f"Trapezoidal Rule of f(x)={f}: {Trapezoidal_Rule(a, b, f, x)}", end="\n\n")
    
    print(f"Newton-Cotes of f(x)={f}: {Newton_Cotes(n, f, x, a, b)}", end="\n\n")
    
    print(f"Cavalieri-Simpson of f(x)={f}: {Cavalieri_Simpson(a, b, f, x)}", end="\n\n")
    
    print(f"Iterated Cavalieri-Simpson of f(x)={f}: {Iterated_Cavalieri_Simpson(a, b, e, n, f, x)[1]}", end="\n\n")
    
    print(f"Iterated Cavalieri-Simpson Error of f(x)={f}: {Iterated_Cavalieri_Simpson_Error(a, b, Iterated_Cavalieri_Simpson(a, b, e, n, f, x)[0], f, x)}", end="\n\n\n")
    
    
    Tn = Itereted_Trapezoidal_Rule(a, b, 1, e, f, x)
    
    print("#"*5, end="")
    print(" Itereted Trapezoidal Rule ", end="")
    print("#"*5)
    
    #print(1/3)
    print(f"Q: {Tn[0]}")
    print(f"Estimated Error: {Tn[1]}")
    print(f"Real Error: {abs((1/3) - Tn[0])}")
    print(f"Points: {Tn[2]}", end="\n\n\n")
    
    g = x**3
    err = float64(10**(-8))
    a, b = float64(0), float64(2*np.pi)
    points = 1
    w = float64(10)
    
    h = cos(w*x)
    
    

    print("#"*5, end="")
    print(" Filon ", end="")
    print("#"*5)
    
    
    #filon = Filon(a, b, points, err, float64(10**(-6)), w, g, x)
    #
    #print(f"Q: {filon[0]}")
    #print(f"Estimated Error: {filon[1]}")
    #print(f"Real Error: {abs(1.18435252813066 - filon[0])}")
    #print(f"h: {filon[2]}")
    #print(f"point: {filon[3]}")
    
    
    f = lambdify(x, f, "numpy")
    
    n = 10000000
    print(1/n)
    yn = np.linspace(0, 1, n)
    
    yn = [f(i) for i in yn]
    
    print(Iterated_Cavalieri_Simpson_With_Points(yn, float64(1/n)))
    