function [X0, tabla, output] = Secante(X0, X1, Tol, N_Max, f, Control_E)

% setup
if Control_E==0
    func_error = @(x, x0) abs(x - x0);
else
    func_error = @(x, x0) abs((x - x0) / x);
end

syms x
Contador = 0;
Funcion_Eval_Tabla(1) = double(subs(f, x, X0));
Funcion_Eval_Tabla(2) = double(subs(f, x, X1));
Error_Tabla(1) = Tol + 1;
Error_Tabla(2) = func_error(X1, X0);
X_n_Tabla(1) = X0;
X_n_Tabla(2) = X1;

while Error_Tabla(Contador + 2) > Tol && Contador < N_Max
    % Aplicar la fórmula de la secante
    X_new = X1 - Funcion_Eval_Tabla(Contador + 2) * (X1 - X0) / (Funcion_Eval_Tabla(Contador + 2) - Funcion_Eval_Tabla(Contador + 1));
    
    X_n_Tabla(Contador + 3) = X_new;
    Funcion_Eval_Tabla(Contador + 3) = double(subs(f, x, X_new));
    Error_Tabla(Contador + 3) = func_error(X_new, X1);
    
    % Actualizar valores
    X0 = X1;
    X1 = X_new;
    Contador = Contador + 1;
end

if Funcion_Eval_Tabla(Contador + 2) == 0
    output = sprintf('%f es una aproximación de una raíz de f(x) con una tolerancia = %f \n', X1, Tol);
elseif Error_Tabla(Contador + 2) < Tol
    output = sprintf('%f es una aproximación de una raíz de f(x) con una tolerancia = %f \n', X1, Tol);
else
    output = sprintf('El método fracasó en %d iteraciones \n', N_Max);
end

Variables_tabla = {'x_n', 'F(x) evaluada', 'Error'};
tabla = table(X_n_Tabla', Funcion_Eval_Tabla', Error_Tabla', 'VariableNames', Variables_tabla);

end
