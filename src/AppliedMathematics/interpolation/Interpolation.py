import numpy as np 
from sympy import *
from sympy.core.function import UndefinedFunction
from scipy import *
from sympy.core.add import Add
from sympy.core.mul import Mul
from sympy.core.function import UndefinedFunction
from numpy import float64

from math import factorial


nullsym = Symbol("")
nullfunc = nullsym
nullfunc = nullfunc.diff(nullsym, 2)
nullfloat = float64(None)



def LS(x: float64, xn: list[float64], yn: list[float64]) -> float64:
    
    y = float64(yn[0] + (x-xn[0])*((yn[0]/(xn[0]-xn[1]) + (yn[1]/(xn[1]-xn[0])))))
    
    return y


def Linear_Spline(x: float64, xn: list[float64]=[], yn: list[float64]=[], f: UndefinedFunction=nullfunc, s: Symbol=nullsym, a: float64=nullfloat, b: float64=nullfloat, h: float64=nullfloat) -> float64:
    
    if len(xn) != len(yn):
        
        raise Exception("Sorry, xn and yn must have the same lenght")
    
    if (x > max(xn) or x < min(xn)) and (len(xn) != 0 and len(yn) != 0):
        
        raise Exception("Sorry x is not in the range")
    
    
    if (f == nullfunc and len(xn) == 0) or (f != nullfunc and (s == nullsym)) or (f != nullfunc and (a == nullfloat or b == nullfloat)):
        
        raise Exception("Not enough data")
    
    if f != nullfunc and len(xn) != 0:
        
        raise Exception("Too much data")
    
    
    if f != nullfunc:
        
        f = lambdify(s, f, "numpy")
        
        n = int((b-a)/h)
        
        xn = np.linspace(a, b, n)
        
        yn = [f(i) for i in xn]
    
    
    
    idx = np.searchsorted(xn, x) - 1
    
    if idx == -1:
        
        idx = len(xn) - 2
    
    return LS(x, xn[idx:idx+2], yn[idx:idx+2])



def Quadratic_Spline(x: float64, xn: list[float64]=[], yn: list[float64]=[], yin: list[float64]=[], f: UndefinedFunction=nullfunc, s: Symbol=nullsym, a: float64=nullfloat, b: float64=nullfloat, h: float64=nullfloat) -> tuple[float64, list[float64]]:
    
    
    if len(xn) != len(yn):
        
        raise Exception("Sorry, xn and yn must have the same lenght")
    
    if (x > max(xn) or x < min(xn)) and (len(xn) != 0 and len(yn) != 0):
        
        raise Exception("Sorry x is not in the range")
    
    
    if (f == nullfunc and len(xn) == 0) or (f != nullfunc and (s == nullsym)) or (f != nullfunc and (a == nullfloat or b == nullfloat)):
        
        raise Exception("Not enough data")
    
    if f != nullfunc and len(xn) != 0:
        
        raise Exception("Too much data")
    
    
    if f != nullfunc:
        
        f = lambdify(s, f, "numpy")
        
        n = int((b-a)/h)
        
        xn = np.linspace(a, b, n)
        
        yn = [f(i) for i in xn]
        
        
    
    if len(xn) < 3:
        
        raise Exception("Not possible")
    
    
    idx = np.searchsorted(xn, x) - 1
    
    
    if idx == 0:
        
        ii = [idx, idx+1, idx+2]
        
        flag = -1
        
    elif idx == len(xn)-1:
        
        ii = [idx-2, idx-1, idx]
        
        flag = 0
        
    else:
        
        ii = [idx-1, idx, idx+1]
        
        flag = 1
        
        
        
    if yin[ii[-1]] == None and yin[ii[0]] == None:
        
        raise Exception("Not enough conditions")
    
    
    
    
    A = [[xn[i]**j for j in range(len(ii))] for i in ii]
    
    
    
    
    for i in range(len(ii)):
        
        A[0].append(0)
        
    for i in range(len(ii)):
        
        A[1].append(0)
        
        
    a = [A[1][i] for i in range(len(ii))]
    
    for i in range(len(ii)):
        
        a.insert(0, 0)
    
    A.insert(2, a)
        
    for i in range(len(ii)):
        
        A[-1].insert(0, 0)
    
    
    
    m1 = []
    
    for j in range(len(ii)):
            
        if xn[ii[1]] == 0 and j-1 < 0:
                
            m1.append(0)
            
        else:
            
            m1.append(j*(xn[ii[1]]**(j-1)))
    
        
    for i in range(len(ii)):
        
        m1.append(m1[i] * -1)
    
    
    A.append(m1)
    
    
    if yin[ii[-1]] != None:
        
    
        A.append([0, 0, 0, 0, 1, 2*xn[ii[-1]]])
        
        b = [yn[idx], yn[ii[1]], yn[ii[1]], yn[ii[-1]], 0, yin[ii[-1]]]
        
        
    else:
        
        A.append([0, 1, 2*xn[ii[0]], 0, 0, 0])
        
        b = [yn[idx], yn[ii[1]], yn[ii[1]], yn[ii[-1]], 0, yin[ii[0]]]
    
    
    an = np.linalg.solve(A, b)
    
    
    qs = 0
    
    
    for i in range(3):
        
        qs += an[i] * (x**i)
    
    
    
     
    return (qs, an)
            
            
            
def Cubic_Spline(x: float64, xn: list[float64]=[], yn: list[float64]=[], yin: list[float64]=[], yiin: list[float64]=[], f: UndefinedFunction=nullfunc, s: Symbol=nullsym, a: float64=nullfloat, b: float64=nullfloat, h: float64=nullfloat) -> tuple[float64, list[float64]]:
    
    if len(xn) != len(yn):
        
        raise Exception("Sorry, xn and yn must have the same lenght")
    
    
    if (x > max(xn) or x < min(xn)) and (len(xn) != 0 and len(yn) != 0):
        
        raise Exception("Sorry x is not in the range")
    
    
    if (f == nullfunc and len(xn) == 0) or (f != nullfunc and (s == nullsym)) or (f != nullfunc and (a == nullfloat or b == nullfloat)):
        
        raise Exception("Not enough data")
    
    if f != nullfunc and len(xn) != 0:
        
        raise Exception("Too much data")
    
    
    if f != nullfunc:
        
        f = lambdify(s, f, "numpy")
        
        n = int((b-a)/h)
        
        xn = np.linspace(a, b, n)
        
        yn = [f(i) for i in xn]
        
        
    
    if len(xn) < 3:
        
        raise Exception("Not possible")
    
    
    idx = np.searchsorted(xn, x) - 1
    
    
    if idx == 0:
        
        ii = [idx, idx+1, idx+2]
        
        flag = -1
        
    elif idx == len(xn)-1:
        
        ii = [idx-2, idx-1, idx]
        
        flag = 0
        
    else:
        
        ii = [idx-1, idx, idx+1]
        
        flag = 1
        
    
    
    A = [[xn[i]**j for j in range(len(ii)+1)] for i in ii]
            
    
    
    for i in range(len(ii)+1):
        
        A[0].append(0)
        
    for i in range(len(ii)+1):
        
        A[1].append(0)
        
        
    a = [A[1][i] for i in range(len(ii)+1)]
    
    for i in range(len(ii)+1):
        
        a.insert(0, 0)
    
    A.insert(2, a)
        
    for i in range(len(ii)+1):
        
        A[-1].insert(0, 0)
        
        
    
    m1 = []

        
    for j in range(len(ii)+1):
            
        if xn[ii[1]] == 0 and j-1 < 0:
                
            m1.append(0)
            
        else:
            
            m1.append(j*(xn[ii[1]]**(j-1)))
    
        
    for i in range(len(ii)+1):
        
        m1.append(m1[i] * -1)

    
    A.append(m1)
        
    m = [[0, 0, 2, 6*xn[i]] for i in ii]
    
    
    for i in range(len(ii)+1):
        
        m[0].append(0)
        
    for i in range(len(ii)+1):
        
        m[1].append(m[1][i] * -1)
        
    for i in range(len(ii)+1):
        
        m[2].insert(0, 0)
        
    
    for i in m:
        
        A.append(i)
        
        
    b = [yn[ii[0]], yn[ii[1]], yn[ii[1]], yn[ii[2]], 0, yiin[ii[0]], 0, yiin[ii[2]]]
    
    
    an = np.linalg.solve(A, b)
    
    
    if flag == -1 or flag == 1:
        
        cs = 0.0
        
        for i in range(0, 4):
            
            cs += an[i]*(x**i)
        
    else:

        cs = 0.0
        
        for i in range(4, len(an)):
            
            cs += an[i]*(x**i)
            
            
            
    return (cs, an)




def Linear_Spline_Error(a: float, b: float, f: UndefinedFunction, s: Symbol, xn: list[float64]=[], h: float=0) -> float64:    
    
    if (len(xn) == 0 and h <= 0) or (f == nullfunc) or (s == nullsym):
        
        raise Exception("Not enough data")
    
    if len(xn) != 0 and h > 0:  
        
        raise Exception("Too much data") 
        
    
    if len(xn) != 0:
        
        for i in range(len(xn)-1):
            
            if abs(xn[i+1]-xn[i]) > h:
                
                h = abs(xn[i+1]-xn[i])
    
    
    
    D = np.linspace(a, b, 10000)
    
    f = f.diff(s, 2)
    
    f = lambdify(s, f, "numpy")
    
    D = max([abs(f(i)) for i in D])
    
    return (D*(h**2))/float64(8)
    
    

def Quadratic_Spline_Error(a: float, b: float, f: UndefinedFunction, s: Symbol, xn: list[float64]=[], h: float=0) -> float64:
    
    if len(xn) == 0 and h <= 0:
        
        raise Exception("Not enough data")
    
    if len(xn) != 0 and h > 0:  
        
        raise Exception("Too much data") 
        
    
    if len(xn) != 0:
        
        for i in range(len(xn)-1):
            
            if abs(xn[i+1]-xn[i]) > h:
                
                h = abs(xn[i+1]-xn[i])
                
                
    D = np.linspace(a, b, 10000)
    
    f = f.diff(s, 3)
    
    f = lambdify(s, f, "numpy")
    
    D = max([abs(f(i)) for i in D])
    
    C = float64(1)/(float64(72)*np.sqrt(3))
    
    return ((h**3) * C * D)


def Cubic_Spline_Error(a: float, b: float, f: UndefinedFunction, s: Symbol, xn: list[float64]=[], h: float=0, constrained: bool=False) -> float64:
    
    if len(xn) == 0 and h <= 0:
        
        raise Exception("Not enough data")
    
    if len(xn) != 0 and h > 0:  
        
        raise Exception("Too much data") 
        
    
    if len(xn) != 0:
        
        for i in range(len(xn)-1):
            
            if abs(xn[i+1]-xn[i]) > h:
                
                h = abs(xn[i+1]-xn[i])
                
                
    D = np.linspace(a, b, 10000)
    
    f = f.diff(s, 4)
    
    f = lambdify(s, f, "numpy")
    
    D = max([abs(f(i)) for i in D])
    
    if constrained:
        
        C = float64(1) / float64(16)
        
    else:
        
        C = float64(5) / float64(384)
        
    
    return (C * D * (h**4))



def Taylor(f: UndefinedFunction, s: Symbol, x: float64, x0: float64, points: int) -> float64:
    
    n = points
    
    tmp = lambdify(s, f, 'numpy')
    der = f.diff(s)
    d = lambdify(s, der, 'numpy')
    
    t = tmp(x0)
    
    
    for i in range(1, n+1):
        
        t += float64(d(x0) * (((x-x0)**(i))/factorial((i))))
        
        der = der.diff(s)
        d = lambdify(s, der, 'numpy')
    
    return t


def Symbolic_Taylor(f: UndefinedFunction, s: Symbol, x: Symbol, x0: float64, points: int) -> Add:
    
    
    n = points
    
    tmp = lambdify(s, f, 'numpy')
    der = f.diff(s)
    d = lambdify(s, der, 'numpy')
    
    t = tmp(x0)
    
    
    for i in range(1, n+1):
        
        t += d(x0) * (((x-x0)**(i))/factorial((i)))
        
        der = der.diff(s)
        d = lambdify(s, der, 'numpy')
    
    return t



def L(vs: list[float64], v: float64, i: int) -> float64:
    
    l = float64(1)
    
    for j in range(len(vs)):
        
        if j != i:
            
            l *= float64((v-vs[j]) / (vs[i]-vs[j]))
            
    return l


def Lagrange(x: float64, xn: list[float64]=[], yn: list[float64]=[], f: UndefinedFunction=nullfunc, s: Symbol=nullsym, a: float64=nullfloat, b: float64=nullfloat, h: float64=nullfloat) -> float64:
    
    if len(xn) != len(yn):
        
        raise Exception("Sorry, xn and yn must have the same lenght")
    
    
    if (len(xn) != 0 and len(yn) != 0) and (x > max(xn) or x < min(xn)):
        
        raise Exception("Sorry x is not in the range")
    
    
    if (f == nullfunc and len(xn) == 0) or (f != nullfunc and (s == nullsym)) or (f != nullfunc and (a == nullfloat or b == nullfloat)):
        
        raise Exception("Not enough data")
    
    if f != nullfunc and len(xn) != 0:
        
        raise Exception("Too much data")
    
    
    if f != nullfunc:
        
        f = lambdify(s, f, "numpy")
        
        n = int((b-a)/h)
        
        xn = np.linspace(a, b, n)
        
        yn = [f(i) for i in xn]
    
    
    l = float64(0)
    
    n = len(xn)
    
    for i in range(n):
        
        l += float64(yn[i] * L(xn, x, i))
        

    return l


def Symbolic_L(vs: list[float64], s: Symbol, i: int) -> Add:
    
    l = float64(1)
    
    for j in range(len(vs)):
        
        if j != i:
            
            l *= ((s-vs[j]) / (vs[i]-vs[j]))
            
    return l



def Symbolic_Lagrange(xn: list[float64]=[], yn: list[float64]=[], f: UndefinedFunction=nullfunc, s: Symbol=nullsym, a: float64=nullfloat, b: float64=nullfloat, h: float64=nullfloat) -> Add:
    
    
    if len(xn) != len(yn):
        
        raise Exception("Sorry, xn and yn must have the same lenght")
    
    
    if (x > max(xn) or x < min(xn)) and (len(xn) != 0 and len(yn) != 0):
        
        raise Exception("Sorry x is not in the range")
    
    
    if (f == nullfunc and len(xn) == 0) or (f != nullfunc and (s == nullsym)) or (f != nullfunc and (a == nullfloat or b == nullfloat)):
        
        raise Exception("Not enough data")
    
    if f != nullfunc and len(xn) != 0:
        
        raise Exception("Too much data")
    
    
    if f != nullfunc:
        
        f = lambdify(s, f, "numpy")
        
        n = int((b-a)/h)
        
        xn = np.linspace(a, b, n)
        
        yn = [f(i) for i in xn]
    
    
    l = float64(0)
    
    n = len(xn)
    
    for i in range(n):
        
        l += yn[i] * Symbolic_L(xn, s, i)
        

    return l



def Divided_Differences(xn: list[float64]=[], yn: list[float64]=[], f: UndefinedFunction=nullfunc, s: Symbol=nullsym, a: float64=nullfloat, b: float64=nullfloat, h: float64=nullfloat) -> float64:
    
    if len(xn) != len(yn):
        
        raise Exception("Sorry, xn and yn must have the same lenght")
    
    
    if (f == nullfunc and len(xn) == 0) or (f != nullfunc and (s == nullsym)) or (f != nullfunc and (a == nullfloat or b == nullfloat)):
        
        raise Exception("Not enough data")
    
    if f != nullfunc:
        
        raise Exception("Too much data")
    
    
    if f != nullfunc and len(xn) == 0:
        
        f = lambdify(s, f, "numpy")
        
        n = int((b-a)/h)
        
        xn = np.linspace(a, b, n)
        
        yn = [f(i) for i in xn]
        
    
    df = float64(0)
    
    n = len(xn)
    
    for i in range(n):
        
        m = 1
        
        for j in range(n):
            
            if i != j:
                
                m *= float64(xn[i] - xn[j])
                
        df += float64(yn[i] / m)
    
    return df



def Symbolic_Newtonian_Polynomials(xn: list[float64]=[], yn: list[float64]=[], f: UndefinedFunction=nullfunc, s: Symbol=nullsym, a: float64=nullfloat, b: float64=nullfloat, h: float64=nullfloat) -> Add:
    
    if len(xn) != len(yn):
        
        raise Exception("Sorry, xn and yn must have the same lenght")
    
    
    if (x > max(xn) or x < min(xn)) and (len(xn) != 0 and len(yn) != 0):
        
        raise Exception("Sorry x is not in the range")
    
    
    if (f == nullfunc and len(xn) == 0) or (f != nullfunc and (s == nullsym)) or (f != nullfunc and (a == nullfloat or b == nullfloat)):
        
        raise Exception("Not enough data")
    
    if f != nullfunc and len(xn) != 0:
        
        raise Exception("Too much data")
    
    
    if f != nullfunc:
        
        f = lambdify(s, f, "numpy")
        
        n = int((b-a)/h)
        
        xn = np.linspace(a, b, n)
        
        yn = [f(i) for i in xn]
    
    n = len(xn)
    
    p = Divided_Differences(xn=[xn[0]], yn=[yn[0]], f=f, s=s, a=a, b=b, h=h)
    
    for i in range(1, n):
        
        tmp = 1
        
        for j in range(i):
            
            tmp *= (s-xn[j])
            
        p += tmp*Divided_Differences(xn=xn[:i+1], yn=yn[:i+1], f=f, s=s, a=a, b=b, h=h)
        
    return p


def Newtonian_Polynomials(x: float64, xn: list[float64]=[], yn: list[float64]=[], f: UndefinedFunction=nullfunc, s: Symbol=nullsym, a: float64=nullfloat, b: float64=nullfloat, h: float64=nullfloat) -> float64:
    
    if len(xn) != len(yn):
        
        raise Exception("Sorry, xn and yn must have the same lenght")
    
    
    if (x > max(xn) or x < min(xn)) and (len(xn) != 0 and len(yn) != 0):
        
        raise Exception("Sorry x is not in the range")
    
    
    if (f == nullfunc and len(xn) == 0) or (f != nullfunc and (s == nullsym)) or (f != nullfunc and (a == nullfloat or b == nullfloat)):
        
        raise Exception("Not enough data")
    
    if f != nullfunc and len(xn) != 0:
        
        raise Exception("Too much data")
    
    
    if f != nullfunc:
        
        f = lambdify(s, f, "numpy")
        
        n = int((b-a)/h)
        
        xn = np.linspace(a, b, n)
        
        yn = [f(i) for i in xn]
        
    
    n = len(xn)
    
    p = Divided_Differences([xn[0]], [yn[0]])
    
    for i in range(1, n):
        
        tmp = 1
        
        for j in range(i):
            
            tmp *= float64(x-xn[j])
            
        p += float64(tmp*Divided_Differences(xn[:i+1], yn[:i+1]))
        
    return p


def Vandermonde_Det(V: list[list[float64]]) -> float64:
    
    n = len(V)
    
    for i in range(n):
        
        for j in V:
            
            if n != len(j):
                
                raise Exception("Non-Square Matrix")
            
    d = float64(1)
    
    for i in range(n-1):
        
        d2 = float64(1)
        
        for j in range(i+1, n):
            
            d2 *= (V[j][1] - V[i][1])
            
        d *= d2
        
    
    return d
    




if __name__ == "__main__":
    
    y = Symbol("y")
    z = Symbol("z")
    
    f = y
    
    x = float64(3.0)
    
    xn = [0.0, 1, 2, 4, 5]
    yn = [0.0, 1, 2, 4, 5]
    yin = [0, 5, 6, 9, 10]
    yin2 = [None, 0, None, None, None]
    yiin = [22, 3, 7, 8, 1]
    
    
    
    print("#"*10, end="")
    print(" Info ", end="")
    print("#"*10)
    
    print(f"Symbols: {y}, {z}")
    print(f"Function: f(y)={f}")
    print(f"Number to find: {x}")
    print(f"xn: {xn}")
    print(f"yn: {yn}")
    print(f"yin: {yin}")
    print(f"yin2: {yin2}")
    print(f"yiin: {yiin}", end="\n\n\n")
    
    
    print("#"*10, end="")
    print(" Linear Spline", end="")
    print("#"*10)
    
    print(f"Lineal Spline: {Linear_Spline(x, xn, yn)}", end="\n\n\n")
    
    
    print("#"*10, end="")
    print(" Taylor and Symbolic Taylor ", end="")
    print("#"*10)
    
    print(f"Taylor: {Taylor(f, y, x, xn[0], len(xn))}")
    
    t = Symbolic_Taylor(f, y, z, xn[0], len(xn))
    
    t = lambdify(z, t, "numpy")
    
    print(f"Symbolic Taylor: {t(x)}", end="\n\n\n")
    
    
    print("#"*10, end="")
    print(" Lagrange and Symbolic Lagrange ", end="")
    print("#"*10)
    
    
    print(f"Lagrange: {Lagrange(x, xn, yn)}")
    
    l = Symbolic_Lagrange(s=y, xn=xn, yn=yn)
    
    l = lambdify(y, l, "numpy")
    
    print(f"Symbolic Lagrange: {l(x)}", end="\n\n\n")
    
    
    print("#"*10, end="")
    print(" Newtonian and Symbolic Newtonian ", end="")
    print("#"*10)
    
    np = Symbolic_Newtonian_Polynomials(xn=xn, yn=xn, s=y)
    
    
    np = lambdify(y, np, "numpy")
    
    print(f"Newtonian: {Newtonian_Polynomials(x, xn, yn)}")
    
    print(f"Symbolic Newtonian: {np(x)}", end="\n\n\n")
    
    
    print("#"*10, end="")
    print(" Quadratic Spline ", end="")
    print("#"*10)
    
    import numpy as np 
    
    qb = Quadratic_Spline(x=x, xn=xn, yn=yn, yin=yin2)
    
    print(f"Quadratic Spline: {qb[0]}, {qb[1]}", end="\n\n\n")
    
    
    print("#"*10, end="")
    print(" Cubic Spline ", end="")
    print("#"*10)
    
    cb = Cubic_Spline(x, xn, yn, yin, yiin)
    
    print(f"Cubic Spline: {cb[0]}, {cb[1]}", end="\n\n\n")
    
    
    