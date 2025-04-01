function [Resultado, tabla, output] = Punto_Fijo(g, f, X0, Tol, N_Max, Control_E)
    % g: función de punto fijo (expresión simbólica en x)
    % f: función original para evaluar f(x) en cada iteración (expresión simbólica en x)
    % X0: aproximación inicial
    % Tol: tolerancia
    % N_Max: número máximo de iteraciones
    % Control_E: 0 para error absoluto, 1 para error relativo

    syms x
    % Inicialización: almacenar la primera aproximación y el valor de f en esa aproximación
    X_n_Tabla = X0;           % Aproximaciones sucesivas
    Error_Tabla = [1];         % Se llenará con los errores; la primera iteración se marcará con NaN
    f_Tabla = double(subs(f, x, X0));  % Valor de f(x) en la aproximación inicial
    
    Contador = 0;
    currentX = X0;
    
    while Contador < N_Max
        % Evaluar la función g en el punto actual
        X1 = double(subs(g, x, currentX));

        % Calcular el error:
        % En la primera iteración, no hay valor anterior para comparar, por lo que se asigna NaN.

        if Control_E == 0
            Error = abs(X1 - currentX);         % Error absoluto
        else
            Error = abs((X1 - currentX) / X1);    % Error relativo
        end


        % Almacenar la nueva aproximación, el error y el valor de f(x) en la tabla
        X_n_Tabla = [X_n_Tabla; X1];
        Error_Tabla = [Error_Tabla; Error];
        f_Tabla = [f_Tabla; double(subs(f, x, X1))];
        
        % Criterio de parada (para iteraciones posteriores a la primera)
        if Contador > 0 && Error < Tol
            break;
        end
        
        % Actualizar currentX y el contador para la siguiente iteración
        currentX = X1;
        Contador = Contador + 1;
    end
    
    % Resultados finales
    Resultado = X1;
    tabla = table(X_n_Tabla, Error_Tabla, f_Tabla, 'VariableNames', {'x_n', 'Error', 'f(x_n)'});
    output = sprintf('%f es una aproximación con tolerancia %f', X1, Tol);
end
