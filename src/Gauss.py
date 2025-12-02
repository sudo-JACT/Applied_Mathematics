import numpy as np 
from sympy import *
from sympy.core.function import UndefinedFunction
from scipy import *



def Golub_Welsch(nodes: int) -> list[list[float]]:
    
    
    return [[0.0,0.0,0.0,0.0], [0.0,0.0,0.0,0.0,0.0]]




def Finde_Xk_and_Wk(nodes: int) -> list[list[float]]:
    
    
    m = []
    
    
    try:
        
        with open("computed_nodes.csv", "r") as file:
            
            l = file.readlines()
            
            l = [i.rstrip().split(",") for i in l]
    
    
            c = 0
            flag = False
            
            for i in l:
                
                if int(i[0]) == nodes:
                    
                    c += 1
                    
                    k = [float(j) for j in i]
                    
                    k.pop(0)
                    
                    m.append(k)
                    
                if c == 2:
                    
                    flag = True
                    
                    break
        
            
            
            if not(flag):
                
                m.clear()
                
                file.close()
                
                with open("computed_nodes.csv", "a") as file:
                    
                    m = Golub_Welsch(nodes)
                    
                    n = "" + str(nodes) + ","
                    w = "" + str(nodes) + ","
                    
                    for i in m[0]:
                        
                        n += str(i) + ","
                    
                    n = n[:-1]
                    n += "\n"
                        
                        
                        
                    for i in m[1]:
                        
                        w += str(i) + ","
                        
                    w = w[:-1]
                    w += "\n"
                    
                    
                    file.write(n)
                    file.write(w)
                        


            
    except:
        
        with open("computed_nodes.csv", "w") as file:
                    
            m = Golub_Welsch(nodes)
            
            n = "" + str(nodes) + ","
            w = "" + str(nodes) + ","
            
            for i in m[0]:
                
                n += str(i) + ","
            
            n = n[:-1]
            n += "\n"
                
                
                
            for i in m[1]:
                
                w += str(i) + ","
                
            w = w[:-1]
            w += "\n"
            
            
            file.write(n)
            file.write(w)
    
    
        
    return m
            
            




def Gauss_Legendre(nodes: int, f: type[UndefinedFunction], s: Symbol) -> float:
    
    m = Finde_Xk_and_Wk(nodes)
    
    G = 0
    
    n = len(m[0])
    
    g = lambdify(s, f, "numpy")
    
    
    for i in range(n):
        
        G += (m[0][i] * g(m[1][i]))
        
    
    return G



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



def Gauss_Chebyshev_2(n: int, f: type[UndefinedFunction], s: Symbol) -> float:
    
    pass




if __name__ == "__main__":
    
    x = Symbol("x")
    
    f = 1
    
    n = 2
    
    print(Gauss_Chebyshev_1(n, f, x))