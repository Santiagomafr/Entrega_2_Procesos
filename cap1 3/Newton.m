function [X0, tabla, output] = Newton(X0, Tol, N_Max, f, Control_E)
    % Newton: Método de Newton para hallar una raíz de f(x)=0.
    % Entradas:
    %   X0       : Aproximación inicial.
    %   Tol      : Tolerancia para el error.
    %   N_Max    : Número máximo de iteraciones.
    %   f        : Función f(x) (expresión simbólica).
    %   Control_E: 0 para error absoluto, 1 para error relativo.
    %
    % Salidas:
    %   X0     : Aproximación final.
    %   tabla  : Tabla que muestra las aproximaciones, los valores de f(x),
    %            los valores de f'(x) y el error en cada iteración.
    %   output : Mensaje de salida.

    % Configuración del cálculo del error
    if Control_E == 0
        func_error = @(x, x0) abs(x - x0);
    else
        func_error = @(x, x0) abs((x - x0) / x);
    end

    syms x
    Derivada = diff(f);  % Cálculo de la derivada de f(x)
    
    % Inicialización
    Contador = 0;
    
    % Evaluación inicial de f(x) y f'(x) en X0
    X_n_Tabla(Contador+1) = X0;
    Funcion_Eval_Tabla(Contador+1) = double(subs(f, x, X0));
    Derivada_Eval_Tabla(Contador+1) = double(subs(Derivada, x, X0));
    Error_Tabla(Contador+1) = Tol + 1;  % Para garantizar que se ingrese al bucle
    Error = Error_Tabla(Contador+1);

    % Bucle iterativo de Newton
    while Error > Tol && Contador < N_Max
        % Calcular la siguiente aproximación usando la fórmula de Newton:
        % x_{n+1} = x_n - f(x_n)/f'(x_n)
        X_n_Tabla(Contador+2) = X0 - Funcion_Eval_Tabla(Contador+1) / Derivada_Eval_Tabla(Contador+1);
        
        % Evaluar f y f' en la nueva aproximación
        Funcion_Eval_Tabla(Contador+2) = double(subs(f, x, X_n_Tabla(Contador+2)));
        Derivada_Eval_Tabla(Contador+2) = double(subs(Derivada, x, X_n_Tabla(Contador+2)));
        
        % Calcular el error según el tipo de control
        Error_Tabla(Contador+2) = func_error(X_n_Tabla(Contador+2), X0);
        Error = Error_Tabla(Contador+2);
        
        % Actualizar X0 para la siguiente iteración
        X0 = X_n_Tabla(Contador+2);
        Contador = Contador + 1;
    end

    % Generar mensaje de salida según el resultado final
    if Funcion_Eval_Tabla(end) == 0
        output = sprintf('%f es una aproximación de una raíz de f(x) con tolerancia %f.', X0, Tol);
    elseif Error < Tol
        output = sprintf('%f es una aproximación de una raíz de f(x) con tolerancia %f.', X0, Tol);
    elseif Derivada_Eval_Tabla(end) == 0
        output = sprintf('%f es una posible raíz múltiple de f(x).', X0);
    else
        output = sprintf('El método fracasó en %d iteraciones.', N_Max);
    end

    % Crear tabla de resultados con columnas:
    % x_n: aproximación
    % F(x) evaluada: valor de f(x)
    % f''(x) evaluada: valor de la derivada f'(x)
    % Error: error en cada iteración
    Variables_tabla = {'x_n', 'F(x)_evaluada', 'f''(x)_evaluada', 'Error'};
    tabla = table(X_n_Tabla', Funcion_Eval_Tabla', Derivada_Eval_Tabla', Error_Tabla', ...
                  'VariableNames', Variables_tabla);
end
