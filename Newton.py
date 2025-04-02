import numpy as np

def newton_monod_complete():
    """Newton's method for solving Monod kinetics equation with detailed iteration output."""
    # Monod model parameters
    mu_max = 0.5        # Maximum growth rate [h⁻¹]
    K_S = 0.2           # Saturation constant [g/L]
    desired_mu = 0.3    # Target growth rate [h⁻¹]
    
    # Specified initial points
    initial_points = [0.1, 2.0]  # [S0, S1] in g/L
    
    # Method configuration
    tol = 1e-6
    max_iter = 50
    
    # Protected Monod function against negative values
    def monod(S):
        return mu_max * max(S, np.finfo(float).eps) / (K_S + max(S, np.finfo(float).eps))
    
    def f(S):
        return monod(S) - desired_mu
    
    def df(S):
        return mu_max * K_S / (K_S + max(S, np.finfo(float).eps))**2
    
    # Process both initial points
    results = []
    for i, S in enumerate(initial_points, 1):
        print(f"\n=== Run with S{i} = {S:.1f} g/L ===")
        print("\nIter\tS [g/L]\t\tμ(S)\t\tf(S)\t\tError")
        print("="*60)
        
        for iter in range(max_iter + 1):
            # Evaluate function and derivative
            f_val = f(S)
            df_val = df(S)
            mu = monod(S)
            
            # Display current iteration
            if iter == 0:
                print(f"{iter:3d}\t{S:9.6f}\t{mu:9.6f}\t{f_val:9.6f}\t{'--':>10}")
            else:
                print(f"{iter:3d}\t{S:9.6f}\t{mu:9.6f}\t{f_val:9.6f}\t{error:10.6f}")
            
            # Check convergence
            if abs(f_val) < tol:
                print(f"\nConvergence reached in {iter} iterations")
                print(f"Final solution: S = {S:.8f} g/L")
                print(f"μ(S) = {mu:.8f} h⁻¹")
                results.append(S)
                break
            
            # Adaptive step control
            step = f_val / (df_val + np.finfo(float).eps)  # Avoid division by zero
            alpha = 1.0
            
            # Backtracking line search
            for _ in range(5):
                S_temp = S - alpha * step
                if S_temp > 0 and abs(f(S_temp)) < abs(f_val):
                    break
                alpha /= 2
            
            # Update
            error = abs(alpha * step)
            S = S - alpha * step
            
            if iter == max_iter:
                print("\nWarning: Maximum iterations reached")
                results.append(S)
    
    return results[-1] if results else None

# Execute the function
S_steady = newton_monod_complete()