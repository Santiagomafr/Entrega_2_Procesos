function [X0, tabla, output] = Newton_Multiple(X0, Tol, N_Max, f, Control_E)

% Configuración del control de error
if Control_E == 0
    func_error = @(x, x0) abs(x - x0);  % Error absoluto
else
    func_error = @(x, x0) abs((x - x0) / x);  % Error relativo
end

syms x
Derivada_1 = diff(f);  % Primera derivada
Derivada_2 = diff(Derivada_1);  % Segunda derivada

Contador = 0;
Funcion_Eval_Tabla(Contador + 1) = double(subs(f, x, X0));
Funcion_Eval = Funcion_Eval_Tabla(Contador + 1);
Derivada_1_Eval_Tabla(Contador + 1) = double(subs(Derivada_1, x, X0));
Derivada_1_Eval = Derivada_1_Eval_Tabla(Contador + 1);
Derivada_2_Eval_Tabla(Contador + 1) = double(subs(Derivada_2, x, X0));
Derivada_2_Eval = Derivada_2_Eval_Tabla(Contador + 1);

Error_Tabla(Contador + 1) = Tol + 1;
Error = Error_Tabla(Contador + 1);
X_n_Tabla(Contador + 1) = X0;

% Iteración del método de Newton para raíces múltiples
while Error > Tol && Contador < N_Max
    X_n_Tabla(Contador + 2) = X0 - (Funcion_Eval * Derivada_1_Eval) / ...
        ((Derivada_1_Eval)^2 - Funcion_Eval * Derivada_2_Eval);
    
    Funcion_Eval_Tabla(Contador + 2) = double(subs(f, x, X_n_Tabla(Contador + 2)));
    Funcion_Eval = Funcion_Eval_Tabla(Contador + 2);
    
    Derivada_1_Eval_Tabla(Contador + 2) = double(subs(Derivada_1, x, X_n_Tabla(Contador + 2)));
    Derivada_1_Eval = Derivada_1_Eval_Tabla(Contador + 2);
    
    Derivada_2_Eval_Tabla(Contador + 2) = double(subs(Derivada_2, x, X_n_Tabla(Contador + 2)));
    Derivada_2_Eval = Derivada_2_Eval_Tabla(Contador + 2);
    
    Error_Tabla(Contador + 2) = func_error(X_n_Tabla(Contador + 2), X0);
    Error = Error_Tabla(Contador + 2);
    
    X0 = X_n_Tabla(Contador + 2);
    Contador = Contador + 1;
end

% Salida del resultado
if Funcion_Eval == 0
    output = sprintf('%f es una aproximación de una raíz de f(x) con una tolerancia = %f \n', X0, Tol);
elseif Error < Tol
    output = sprintf('%f es una aproximación de una raíz de f(x) con una tolerancia = %f \n', X0, Tol);
elseif Derivada_1_Eval == 0
    output = sprintf('%f es una posible raíz múltiple de f(x) \n', X0);
else
    output = sprintf('Fracasó en %f iteraciones \n', N_Max);
end

% Creación de la tabla de iteraciones
Variables_tabla = {'x_n', 'F(x) evaluada', 'Error'};
tabla = table(X_n_Tabla', Funcion_Eval_Tabla', Error_Tabla', 'VariableNames', Variables_tabla);

end