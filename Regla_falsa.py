import numpy as np

# Monod function 
def monod(S, mu_max=1, K_s=1.5, mu_deseado=0.7):
    return mu_max * S / (K_s + S) - mu_deseado

def false_position(a, b, tol=1e-7, max_iter=50, mu_max=1, K_s=1.5, mu_d=0.7):
    print("\nFalse position method")
    print("Iter |     a     |     b     |     c     |   f(c)    |  Error  ")
    print("-"*70)
    
    for i in range(max_iter):
        fa = monod(a, mu_max, K_s, mu_d)
        fb = monod(b, mu_max, K_s, mu_d)
        c = (a*fb - b*fa)/(fb - fa)
        fc = monod(c, mu_max, K_s, mu_d)
        
        error = abs(b- a)

        
        print(f"{i:3d} | {a:9.6f} | {b:9.6f} | {c:9.6f} | {fc:9.6f} | {error:9.6f}")
        
        if abs(error) < tol:
            print("\n¡Convergence achived!")
            return c
        
        if fa*fc < 0:
            b = c
        else:
            a = c
    
    print("\nMaximum number of iterations achieved ")
    return c
# Parameters
mu_max = 1
K_s = 1.5
mu_d = 0.7

# Range 
a, b =3, 4

# solution
sol_fp = false_position(a, b, mu_max = mu_max, K_s= K_s, mu_d= mu_d)
print(f"\nSolution false position: S = {sol_fp:.2f} g/L")
