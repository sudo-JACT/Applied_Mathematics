import numpy as np 
from sympy import *
from sympy.core.function import UndefinedFunction
from scipy import *



def write_on_file():
    
    pass



def Golub_Welsch(nodes: int) -> list[list[float]]:
    
    
    return [[0.0,0.0,0.0,0.0], [0.0,0.0,0.0,0.0,0.0]]




def Finde_Xk_and_Wk(nodes: int) -> list[list[float]]:  # funzione che restituisce sia i pesi che i nodi di Legendre
        
    
    m = [
        [[2], [0]],
        [[1, 1], [-0.5773502691896258, 0.5773502691896258]],
        [[0.5555555555555556, 0.8888888888888888, 0.5555555555555556], [-0.7745966692414834, 0, 0.7745966692414834]],
        [[0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538], [-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526]],
        [[0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891], [-0.9061798459386640, -0.5384693101056831, 0, 0.5384693101056831, 0.9061798459386640]]
    ]
    
        
    return m[nodes-1]
            
            




def Gauss_Legendre(nodes: int, f: type[UndefinedFunction], s: Symbol, a: float, b: float) -> float:
    
    m = Finde_Xk_and_Wk(nodes)
    
    G = 0
    
    n = len(m[0])
    
    g = lambdify(s, f, "numpy")
    
    
    for i in range(n):
        
        x = (((b-a)/2) * m[1][i]) + ((a+b)/2)
        
        G += (m[0][i] * g(x))
        
    
    return G * ((b-a) / 2)



def Chebyshev_Polynomials_1(degree: int) -> type[UndefinedFunction]:
    
    x = Symbol("x")
    
    
    if degree == 0:
        
        Tn = 1
    
    elif degree == 1:
        
        Tn = x
    
    else:
        
        Tn = ((2 * x * Chebyshev_Polynomials_1(degree-1)) - Chebyshev_Polynomials_1(degree-2))
        
    
    return Tn



def Chebyshev_Zeros_1(n: int, i: int):
    
    xi = np.cos(((((2*i) + 1) * np.pi) / ((2*n))))
    
    return xi




def Gauss_Chebyshev_1(n: int, f: type[UndefinedFunction], s: Symbol) -> float:
    
    C = 0
    
    w = np.pi/(n)
    
    
    g = lambdify(s, f, "numpy")
    
    for i in range(n):
        
        C += g(Chebyshev_Zeros_1(n, i))
        
    return w*C




def Chebyshev_Zeros_2(n: int, i: int):
    
    xi = np.cos(((i+1) * np.pi) / ((n+2)))
    
    return xi


def Chebyshev_W_2(n: int, i: int) -> float:
    
    wi = ((np.pi) / (n+2)) * (np.sin(((i + 1) * np.pi) / (n+2))**2)
    
    return wi


def Gauss_Chebyshev_2(n: int, f: type[UndefinedFunction], s: Symbol) -> float:
    
    C = 0
    
    
    g = lambdify(s, f, "numpy")
    
    
    for i in range(n+1):
        
        C += (g(Chebyshev_Zeros_2(n, i)) * Chebyshev_W_2(n, i))
        
    
    return C




if __name__ == "__main__":
    
    x = Symbol("x")
    
    f = 1

    n = 5
    
    print(Gauss_Chebyshev_2(n, f, x))
    print(Gauss_Legendre(n, f, x, -1, 1))