def L1d(vs: list[float], v: float, i: int) -> float:
    
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
            
            d += c * L1d(ys, y, j) 
        
        pxy += d * L1d(xs, x, i) 
                
    
    return pxy