import numpy as np

def monod(S, mu_max=1, K_s=1.5, mu_deseado=0.7):
    return mu_max * S / (K_s + S) - mu_deseado

def biseccion_con_tabla(a, b, tol=1e-7, max_iter=50, mu_max=1, K_s=1.5, mu_deseado=0.7):

    # Verificar condición inicial
    if monod(a, mu_max, K_s, mu_deseado) * monod(b, mu_max, K_s, mu_deseado) >= 0:
        raise ValueError("El intervalo [a, b] no cumple f(a)*f(b) < 0.")
    
    print("\nMÉTODO DE BISECCIÓN")
    print("Iter |     a     |     b     |     c     |   f(c)    |  Error  ")
    print("-" * 70)
    
    for i in range(max_iter):
        c = (a + b) / 2
        fc = monod(c, mu_max, K_s, mu_deseado)
        error = abs(b - a)
        
        # Imprimir fila de la tabla
        print(f"{i:3d} | {a:9.6f} | {b:9.6f} | {c:9.6f} | {fc:9.6f} | {error:9.6f}")
        
        # Verificar convergencia
        if abs(error) < tol:
            print("\n¡Convergencia alcanzada!")
            return c
        
        # Actualizar intervalo
        if monod(a, mu_max, K_s, mu_deseado) * fc < 0:
            b = c
        else:
            a = c
    
    print("\nAdvertencia: Máximo de iteraciones alcanzado.")
    return c

# Parámetros del modelo
mu_max = 1    # h⁻¹
K_s = 1.5       # g/L
mu_deseado = 0.7  # h⁻¹

a, b = 1, 4 

print(f"Resolviendo μ(S) = {mu_max}·S/({K_s} + S) = {mu_deseado}")

# Ejecutar método de bisección
sol_biseccion = biseccion_con_tabla(a, b, mu_max=mu_max, K_s=K_s, mu_deseado=mu_deseado)
print(f"\nSolución final: S = {sol_biseccion:.6f} g/L")
