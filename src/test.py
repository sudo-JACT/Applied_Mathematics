import numpy as np 
from sympy import *
from sympy.core.function import UndefinedFunction
from scipy import *


from Gauss import Gauss_Chebyshev_1


if __name__ == "__main__":
    
    x = Symbol("x")
    n = 2
    
    f1 = 1
    f2 = x
    f3 = x**2
    f4 = x**3
    f5 = x**4
    
    
    print("Test 1 - f(x) = 1:")
    print(f"Risultato: {Gauss_Chebyshev_1(n, f1, x)} (Atteso: π ≈ {np.pi})")
    
    print("\nTest 2 - f(x) = x:")
    print(f"Risultato: {Gauss_Chebyshev_1(n, f2, x)} (Atteso: 0)")
    
    print("\nTest 3 - f(x) = x^2:")
    print(f"Risultato: {Gauss_Chebyshev_1(n, f3, x)} (Atteso: π/2 ≈ {np.pi/2})")
    
    print("\nTest 4 - f(x) = x^3 (grado 3, ancora esatto con n=2):")
    
    print(f"Risultato: {Gauss_Chebyshev_1(n, f4, x)} (Atteso: 0)")

    print(f"\nTest 5 - f(x) = x^4: {Gauss_Chebyshev_1(n, f5, x)}")
    