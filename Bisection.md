function [result_table, S_root] = bisection(mu_max, K_s, desired_mu, a, b, nmax, tol)
    f = @(S) mu_max * (S / (K_s + S)) - desired_mu;
    
    % Interval validation
    if f(a) * f(b) >= 0
        error('Error: No sign change in [a, b].');
    end
    
    % Initialization
    table_data = zeros(nmax, 4);
    iter = 0;
    S_root = NaN;
    
    % Main loop
    while iter < nmax
        iter = iter + 1;
        S = (a + b) / 2;
        f_S = f(S);
        current_error = abs(b - a);
        
        % Store results
        table_data(iter, :) = [iter, S, f_S, current_error];
        
        % Stopping criterion
        if current_error <= tol
            S_root = S;
            break;
        end
        
        % Update interval
        if f(a) * f_S < 0
            b = S;
        else
            a = S;
        end
    end
    
    % Format output
    table_data = table_data(1:iter, :);
    headers = {'Iteration', 'S', 'f(S)', 'Error'};
    result_table = array2table(table_data, 'VariableNames', headers);
    
    % Warning if no convergence
    if isnan(S_root)
        warning('No convergence in %d iterations. Last error: %f', nmax, current_error);
        S_root = S;
    end
end

