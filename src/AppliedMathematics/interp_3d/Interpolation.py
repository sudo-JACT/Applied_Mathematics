import numpy as np 
from sympy import lambdify, Symbol
from sympy.core.function import UndefinedFunction
from numpy import float64

from ..core import nullfloat, nullfunc, nullsym
from ..interpolation import Divided_Differences


def L1d(vs: list[float64], v: float64, i: int) -> float64:
    
    l = float64(1)
    
    for j in range(len(vs)):
        
        if j != i:
            
            l *= float64((v-vs[j]) / (vs[i]-vs[j]))
            
    return l




def Lagrangian_R2(xs: list[float64], ys: list[float64], x: float64, y: float64, fxy: UndefinedFunction=nullfunc, s: Symbol=nullsym, t: Symbol=nullsym, zs: list[list[float64]]=[]) -> float64 :  # funzione che calcola il polinomio di Lagrange in R2
    
    
    if (fxy == nullfunc and len(zs) == 0) or (fxy != nullfunc and (s == nullsym or t == nullsym)):
        
        raise Exception("Not enough data")
    
    if fxy != nullfunc and len(zs) != 0:
        
        raise Exception("Too much data")
    
    
    
    pxy = float64(0)
    
    fxy = lambdify([s, t], fxy, "numpy")
    
    
    if len(zs) == 0:
    
        for i in range(len(xs)):
            
            d = float64(0)
            
            for j in range(len(ys)):
                
                c = fxy(xs[i], ys[j])
                
                d += c * L1d(ys, y, j) 
            
            pxy += d * L1d(xs, x, i) 
            
    else:
        
        for i in range(len(xs)):
            
            d = float64(0)
            
            for j in range(len(ys)):
                
                c = zs[i][j]
                
                d += c * L1d(ys, y, j) 
            
            pxy += d * L1d(xs, x, i) 
                
    
    return pxy



def DivDif_2(xn: list[float64], yn: list[float64], zn: list[list[float64]], n: int, m: int) -> float64:
    
    d = float64(0)
    
    
    for i in range(n+1):
        
        
        f = float64(0)
        
        
        
        
        for j in range(m+1):
            
            
            d2 = float64(1)
            
            for k in range(m+1):
                
                if k != j:
                    
                    d2 *= (yn[j] - yn[k])
                    
                    
            f += zn[i][j] / d2
            
            
        d3 = float64(1)
        
        for k in range(n+1):
            
            if k != i:
                
                d3 *= (xn[i] - xn[k])
                
        
        d += f/d3
                    
    
    return d



def Weights(x: float64, xn: list[float64], k: int) -> float64:
    
    w = float64(1)
    
    if k == 0:
        
        return w
    
    
    for i in range(k):
    
        w *= (x-xn[i])
        
    return w




def Newtonian_Polynomials_R2(xn: list[float64], yn: list[float64], x: float64, y: float64, zn: list[list[float64]]=[], f: UndefinedFunction=nullfunc, s: Symbol=nullsym, t: Symbol=nullsym) -> float64:
    
    
    if (len(zn) == 0) and (f == nullfloat or (s == nullsym or t == nullsym)):
        
        raise Exception("Not enough data")
    
    
    if len(zn) == 0:
        
        f = lambdify([s, t], f, "numpy")
        
        zn = [[f(i, j) for j in yn] for i in xn]
    
    
    
    
    d = float64(0)
    
    n = len(xn)
    m = len(yn)
    
    for i in range(n):
        
        wx = Weights(x, xn, i)
        
        d2 = float64(0)
        
        for j in range(m):
            
            wy = Weights(y, yn, j) 
            
            d2 += wy * DivDif_2(xn, yn, zn, i, j)
            
            
        d += wx * d2
    
    
    
        
    return d





if __name__ == "__main__":
    
    import matplotlib.pyplot as plt
    
    
    x = Symbol('x')  # definisco x come symbol per utilizzarlo nelle espressioni matematiche 
    y = Symbol('y')
    fxy = (x**3)-(y**2)
    

    xmin = -1  # definisco a, b, c, d
    xmax = 1
    
    ymin = -2
    ymax = 2
    
    n = 20
    
    interpolation_n = int(n/4)
    

    xs = np.linspace(xmin, xmax, n)  # creo le decomposizioni di passo 1/n
    ys = np.linspace(ymin, ymax, n)  
            
        
    zx = np.linspace(xmin, xmax, interpolation_n)
    zy = np.linspace(ymin, ymax, interpolation_n)
    
    zx2 = np.linspace(xmin, xmax, 2*interpolation_n)
    zy2 = np.linspace(ymin, ymax, 2*interpolation_n)
    
    zs = np.array([[Lagrangian_R2(xs=zx, ys=zy, fxy=fxy, s=x, t=y, x=xi, y=yi, zs=[]) for yi in ys] for xi in xs])  # calcolo il polinomio du 5 punti e poi il polinomio su 10
    
    zs2 = np.array([[Lagrangian_R2(xs=zx2, ys=zy2, fxy=fxy, s=x, t=y, x=xi, y=yi, zs=[]) for yi in ys] for xi in xs])
        
    
    
    fxy = lambdify(["x", "y"], fxy, "numpy")
    
    
    X, Y = np.meshgrid(xs, ys, indexing="ij")  # creo le mesh per i diversi grafici
    
    zx, zy = np.meshgrid(zx, zy, indexing="ij")
    
    zx2, zy2 = np.meshgrid(zx2, zy2, indexing="ij")
    
    Z_real = fxy(X, Y)   # calcolo la funzione sulla mesh
    
    
    
    fig = plt.figure(figsize=(10, 8))  # plotto tutti i grafici 
    ax1 = plt.axes(projection='3d')
    
    surf1 = ax1.plot_surface(X, Y, Z_real, cmap='viridis', alpha=0.9, linewidth=0, antialiased=True)
    ax1.set_title('Funzione Reale: $x^3 - y^2$', fontsize=14, fontweight='bold')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_zlabel('z')
    
    
    
    
    fig = plt.figure(figsize=(10, 8))
    ax2 = plt.axes(projection='3d')
    
    surf2 = ax2.plot_surface(X, Y, zs, cmap='magma', alpha=0.9, linewidth=0, antialiased=True)
    ax2.set_title(f'Interpolazione di Lagrange interpolato su {interpolation_n} punti', fontsize=14, fontweight='bold')
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    ax2.set_zlabel('z')
    
    
    
    fig = plt.figure(figsize=(10, 8))
    ax5 = plt.axes(projection='3d')
    
    surf5 = ax5.plot_surface(X, Y, zs2, cmap='Greys', alpha=0.9, linewidth=0, antialiased=True)
    ax5.set_title(f'Interpolazione di Lagrange interpolato su {2*interpolation_n} punti', fontsize=14, fontweight='bold')
    ax5.set_xlabel('x')
    ax5.set_ylabel('y')
    ax5.set_zlabel('z')
    
    
    
    
    fig = plt.figure(figsize=(10, 8))
    ax3 = plt.axes(projection='3d')
    
    
    e1 = fxy(zx, zy)  # calcolo e plotto l'errore per il polinomio interpolato su 5 punti
    error = np.abs((e1) - zs)
    surf3 = ax3.plot_surface(zx, zy, error, cmap='inferno', alpha=0.9, linewidth=0, antialiased=True)
    ax3.set_title(f'Errore Assoluto di Lagrange interpolato su {interpolation_n} punti', fontsize=14, fontweight='bold')
    ax3.set_xlabel('x')
    ax3.set_ylabel('y')
    ax3.set_zlabel('Errore')
    
    
    
    fig = plt.figure(figsize=(10, 8))
    ax4 = plt.axes(projection='3d')
    
    
    e2 = fxy(zx2, zy2)  # calcolo e plotto l'errore per il polinomio interpolato su 10 punti
    error = np.abs((e2) - zs2)
    surf4 = ax4.plot_surface(zx2, zy2, error, cmap='Greys', alpha=0.9, linewidth=0, antialiased=True)
    ax4.set_title(f'Errore Assoluto di Lagrange interpolato su {2*interpolation_n} punti', fontsize=14, fontweight='bold')
    ax4.set_xlabel('x')
    ax4.set_ylabel('y')
    ax4.set_zlabel('Errore')
    
    
    vmin = min(np.min(Z_real), np.min(zs))
    vmax = max(np.max(Z_real), np.max(zs))
    surf1.set_clim(vmin, vmax)
    surf2.set_clim(vmin, vmax)
    surf5.set_clim(vmin, vmax)

    # Aggiungi colorbar
    fig.colorbar(surf1, ax=ax1, shrink=0.6, aspect=20, label='Valore z')
    fig.colorbar(surf2, ax=ax2, shrink=0.6, aspect=20, label='Valore z')
    fig.colorbar(surf5, ax=ax5, shrink=0.6, aspect=20, label='Valore z')
    fig.colorbar(surf3, ax=ax3, shrink=0.6, aspect=20, label='Errore')
    fig.colorbar(surf4, ax=ax4, shrink=0.6, aspect=20, label='Errore')

    plt.tight_layout()
    
    print("\n\n")  # printo in consolo l'errore dei due polinomi interpolatori 
    print(f"Errore Assoluto di Lagrange interpolato su {interpolation_n} punti = {np.max(np.abs((e1) - zs)):.20f}")
    print(f"\nErrore Assoluto di Lagrange interpolato su {2*interpolation_n} punti = {np.max(np.abs((e2) - zs2)):.20f}")
    
    
    plt.show()