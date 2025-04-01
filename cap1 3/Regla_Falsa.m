function [Resultado, tabla, output] = Regla_Falsa(a, b, Tol, N_Max, Funcion_Objetivo, Control_E)
    syms x
    Contador = 0;
    X_n_Tabla = [];
    Error_Tabla = [];
    Funcion_Eval_Tabla = [];

    % Evaluamos la función en los extremos
    fa = double(subs(Funcion_Objetivo, x, a));
    fb = double(subs(Funcion_Objetivo, x, b));

    % Verificamos que el intervalo sea válido
    if fa * fb > 0
        error('La función debe cambiar de signo en el intervalo [a, b].');
    end

    % Iteración del método
    while Contador < N_Max
        % Nueva aproximación de la raíz
        c = (a * fb - b * fa) / (fb - fa);
        fc = double(subs(Funcion_Objetivo, x, c));

        % Guardamos valores en tablas
        X_n_Tabla = [X_n_Tabla; c];
        Funcion_Eval_Tabla = [Funcion_Eval_Tabla; fc];

        % Cálculo del error
        if Contador == 0
            Error = Tol + 1;  % Inicializamos con un valor mayor a Tol
        else
            if Control_E == 0
                Error = abs(c - X_n_Tabla(end - 1));  % Error absoluto
            else
                Error = abs((c - X_n_Tabla(end - 1)) / c);  % Error relativo
            end
        end
        Error_Tabla = [Error_Tabla; Error];

        % Criterios de parada
        if abs(fc) < Tol || Error < Tol
            break;
        end

        % Actualización del intervalo
        if fa * fc < 0
            b = c;
            fb = fc;
        else
            a = c;
            fa = fc;
        end

        Contador = Contador + 1;
    end

    % Resultados
    Resultado = c;
    tabla = table(X_n_Tabla, Funcion_Eval_Tabla, Error_Tabla, 'VariableNames', {'x_n', 'F(x)', 'Error'});
    output = sprintf('%f es una aproximación con tolerancia %f', c, Tol);
end
