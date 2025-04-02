import numpy as np

# Definition of the Monod function 
def f(S, mu_max=1, K_s=1.5, mu_deseado=0.7):
    return mu_max * (S / (K_s + S)) - mu_deseado

def secante_table(S0, S1, tol=1e-7, max_iter=50, mu_max=1, K_s=1.5, mu_deseado=0.7):
    print("\nIteration |     S0     |     S1     |    f(S1)    |   Error  ")
    print("-" * 60)
    
    for i in range(max_iter):
        f_S1 = f(S1, mu_max, K_s, mu_deseado)
        f_S0 = f(S0, mu_max, K_s, mu_deseado)
        
        # This line used to calculate the realtive errror 
        error = abs(S1 - S0)/S1 if i > 0 else np.nan
        
        # This line is used to print the table 
        print(f"{i:5d}    | {S0:10.9f} | {S1:10.9f} | {f_S1:10.9f} | {error:10.7f}")
        
        if abs(f_S1) <= tol:
            print("\n¡Convergence achieved!")
            return S1
        
        denominator = f_S1 - f_S0
        if abs(denominator) <= 1e-12:
            raise ValueError("¡Division by zero! changes the starting points.")
        
        S_next = S1 - f_S1 * (S1 - S0) / denominator
        S0, S1 = S1, S_next
    
    print("\nWarning: maximum number of iterations reached.")
    return S1

# Parameters and execution 
mu_max = 1
K_s = 1.5
mu_deseado = 0.7
S0, S1 = 0.8, 2.0

print(f"Monod model: μ_max = {mu_max}, K_s = {K_s}, μ_deseado = {mu_deseado}\n")
sol_secante = secante_table(S0, S1, mu_max=mu_max, K_s=K_s, mu_deseado=mu_deseado)
print(f"\nFinal Solution: S = {sol_secante:.6f} g/L")
