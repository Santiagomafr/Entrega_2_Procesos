function [X0, tabla, output] = Newton_Multiple(X0, Tol, N_Max, f, Control_E)
    % Newton_Multiple: Método de Newton modificado para raíces múltiples.
    % Entradas:
    %   X0       : Aproximación inicial.
    %   Tol      : Tolerancia para el error.
    %   N_Max    : Número máximo de iteraciones.
    %   f        : Función f(x) (expresión simbólica).
    %   Control_E: 0 para error absoluto, 1 para error relativo.
    %
    % Salidas:
    %   X0     : Aproximación final.
    %   tabla  : Tabla que muestra x_n, f(x_n), f'(x_n), f''(x_n) y el error en cada iteración.
    %   output : Mensaje de salida.
    
    if Control_E == 0
        func_error = @(x, x0) abs(x - x0);  % Error absoluto
    else
        func_error = @(x, x0) abs((x - x0) / x);  % Error relativo
    end

    syms x
    Derivada_1 = diff(f);      % Primera derivada
    Derivada_2 = diff(Derivada_1);  % Segunda derivada

    Contador = 0;
    
    % Evaluaciones iniciales en X0
    X_n_Tabla(Contador + 1) = X0;
    Funcion_Eval_Tabla(Contador + 1) = double(subs(f, x, X0));
    Derivada_1_Eval_Tabla(Contador + 1) = double(subs(Derivada_1, x, X0));
    Derivada_2_Eval_Tabla(Contador + 1) = double(subs(Derivada_2, x, X0));
    Error_Tabla(Contador + 1) = 1;  % No hay error en la primera iteración
    Error = Tol + 1;  % Para iniciar el bucle

    while Error > Tol && Contador < N_Max
        % Fórmula modificada de Newton para raíces múltiples:
        X_n_Tabla(Contador + 2) = X0 - (Funcion_Eval_Tabla(Contador + 1) * Derivada_1_Eval_Tabla(Contador + 1)) / ...
            ((Derivada_1_Eval_Tabla(Contador + 1))^2 - Funcion_Eval_Tabla(Contador + 1) * Derivada_2_Eval_Tabla(Contador + 1));
        
        % Evaluar f, f' y f'' en la nueva aproximación
        Funcion_Eval_Tabla(Contador + 2) = double(subs(f, x, X_n_Tabla(Contador + 2)));
        Derivada_1_Eval_Tabla(Contador + 2) = double(subs(Derivada_1, x, X_n_Tabla(Contador + 2)));
        Derivada_2_Eval_Tabla(Contador + 2) = double(subs(Derivada_2, x, X_n_Tabla(Contador + 2)));
        
        % Calcular el error
        Error_Tabla(Contador + 2) = func_error(X_n_Tabla(Contador + 2), X0);
        Error = Error_Tabla(Contador + 2);
        
        % Actualizar la aproximación y el contador
        X0 = X_n_Tabla(Contador + 2);
        Contador = Contador + 1;
    end

    % Generar mensaje de salida
    if Funcion_Eval_Tabla(end) == 0
        output = sprintf('%f es una aproximación de una raíz de f(x) con tolerancia = %f.', X0, Tol);
    elseif Error < Tol
        output = sprintf('%f es una aproximación de una raíz de f(x) con tolerancia = %f.', X0, Tol);
    elseif Derivada_1_Eval_Tabla(end) == 0
        output = sprintf('%f es una posible raíz múltiple de f(x).', X0);
    else
        output = sprintf('Fracasó en %d iteraciones.', N_Max);
    end

    % Crear la tabla de resultados con las columnas: x_n, f(x) evaluada, f'(x) evaluada, f''(x) evaluada y Error.
    Variables_tabla = {'x_n', 'F(x)_evaluada', 'f''(x)_evaluada', 'f''''(x)_evaluada', 'Error'};
    tabla = table(X_n_Tabla', Funcion_Eval_Tabla', Derivada_1_Eval_Tabla', Derivada_2_Eval_Tabla', Error_Tabla', ...
                  'VariableNames', Variables_tabla);
end
