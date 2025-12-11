from math import factorial
import numpy as np 
from sympy import *
from sympy.core.function import UndefinedFunction



def binomial(k: int, n: int) -> float:
    
    return (factorial(n) / (factorial(k) * factorial(n-k)))




def Bernstein_Polynomials(k: int, n: int, x: float) -> float :
    
    return (binomial(k, n) * (x**k) * ((1-x)**(n-k)))



def Bernstein_Approximation(k: int, n: int, x: float, f, s) -> float :  # funzione che calcola i polinomi di Bernstein Bf(x)
    
    f = lambdify(s, f, "numpy")
    
    bf = 0
    
    for i in range(k, n+1):
    
        bf += f(i/n) * Bernstein_Polynomials(i, n, x)
        
    return bf







def Bezier_Curves(p: list[float], t: float) -> float:  # curva di Bézier
    
    b = 0
    
    n = len(p)

    for k in range(n):
        
        b += p[k] * Bernstein_Polynomials(k, n-1, t)
        
        
    return b





if __name__ == "__main__":
    
    import matplotlib.pyplot as plt
    
    h = 100
    
    t = [i/h for i in range(h+1)]  # decomposizione t che andara da 0 a 1 con passo 0.01
    
    
    xsup = [0.0, 0.004712, 0.007052, 0.012157, 0.024309, 0.049040, 0.073920, 0.098786, 0.148638, 0.198551, 0.248554, 0.298540, 0.398560, 0.498566, 0.598899, 1.0]
    ysup = [0.0, 0.012211, 0.014766, 0.018787, 0.029274, 0.040657, 0.045745, 0.051431, 0.057708, 0.061386, 0.061266, 0.061846, 0.061006, 0.060766, 0.046634, 0.0]
    
    xsup = np.array(xsup)
    ysup = np.array(ysup)
    
    
    Bsupx = [Bezier_Curves(xsup, i) for i in t]  # calcolo Bézier per le x superiori
    Bsupy = [Bezier_Curves(ysup, i) for i in t]  # calcolo Bézier per le y superiori
    
    Bsupx = np.array(Bsupx)
    Bsupy = np.array(Bsupy)
    
    
    #plotto le curve
        
    plt.plot(Bsupx, Bsupy, "k", label='Curva di Bézier superiore')
    
    plt.plot(xsup, ysup, "ok", label='Punti Controllo Sup.')
    
    
    
    xinf = [0.0, 0.005209, 0.007681, 0.012940, 0.025467, 0.050649, 0.075755, 0.100876, 0.151044, 0.201198, 0.251319, 0.301411, 0.401417, 0.501361, 0.601185, 1.0]
    yinf = [0.0, -0.008877, -0.011919, -0.014395, -0.019799, -0.027505, -0.032012, -0.037119, -0.044235, -0.050752, -0.055869, -0.059787, -0.060027, -0.057668, -0.050212, 0.0]
    
    Binfx = [Bezier_Curves(xinf, i) for i in t] # calcolo Bézier per le x inferiori
    Binfy = [Bezier_Curves(yinf, i) for i in t] # calcolo Bézier per le y superiori
    
    
    #plotto sia le curve che i punti
        
    plt.plot(Binfx, Binfy, 'r', label='Curva di Bézier inferiore')
    
    plt.plot(xinf, yinf, 'or', label='Punti Controllo Inf.')
    
    
    plt.axis('equal') 
    plt.title('Profilo Alare Boing 707 (Curva di Bézier)')
    plt.legend()
    plt.grid(True)
    
    plt.figure()
    
    plt.plot(Binfx, Binfy, 'r', label='Curva di Bézier inferiore')
    plt.plot(Bsupx, Bsupy, "k", label='Curva di Bézier superiore')


    plt.axis('equal') 
    plt.title('Profilo Alare Boing 707 (Curva di Bézier)')
    
    plt.legend()
    
    plt.grid(True)
    plt.show()