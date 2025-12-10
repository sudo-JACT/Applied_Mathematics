from math import factorial

def Horner(a: list[float], x: float) -> float:
    
    p = a[-1]
    
    for i in range(len(a)-1):
        
        p = (p * x) + a[((len(a)-2) - i)]
    
    return p



def binomial(k: int, n: int) -> float:
    
    return (factorial(n) / (factorial(k) * factorial(n-k)))


if __name__ == "__main__":
    
    a = [1, 2, 3, 7.0]
    
    x = 2
    
    print(Horner(a, x))