function [Resultado, tabla, output] = Biseccion(a, b, Tol, N_Max, Funcion_Objetivo, Control_E)
    syms x
    Contador = 0;
    X_n_Tabla = [];
    Error_Tabla = [];
    Funcion_Eval_Tabla = [];

    % Verificación de cambio de signo
    fa = double(subs(Funcion_Objetivo, x, a));
    fb = double(subs(Funcion_Objetivo, x, b));

    if fa * fb >= 0
        output = 'El intervalo no cumple el criterio de cambio de signo';
        Resultado = NaN;
        tabla = table();
        return;
    end

    while (b - a) / 2 > Tol && Contador < N_Max
        c = (a + b) / 2;
        fc = double(subs(Funcion_Objetivo, x, c));
        
        % Guardar valores en las tablas
        X_n_Tabla = [X_n_Tabla; c];
        Funcion_Eval_Tabla = [Funcion_Eval_Tabla; fc];
        
        % Cálculo del error según el tipo de control
        if Contador == 0
            Error = Tol + 1;  % Inicialización del error para la primera iteración
        else
            if Control_E == 0
                Error = abs(b - a) / 2;  % Error absoluto
            else
                Error = abs((b - a) / c);  % Error relativo
            end
        end
        Error_Tabla = [Error_Tabla; Error];

        % Criterio de parada
        if fc == 0 || Error < Tol
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
    tabla = table(X_n_Tabla, Funcion_Eval_Tabla,  Error_Tabla, 'VariableNames', {'x_n', 'F(x)', 'Error'});
    output = sprintf('%f es una aproximación de una raíz con tolerancia %f', c, Tol);
end
