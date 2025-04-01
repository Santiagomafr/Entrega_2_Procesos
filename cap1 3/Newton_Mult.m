function [Resultado,tabla,output]=Newton_Mult(f,Punto_Inicial,tolerancia,m,niter,Control_E)
syms x
Derdivada=diff(f);
Puntos=[Punto_Inicial];
Punto_Nuevo=Punto_Inicial-m*(Funcion_eval(Punto_Inicial,f)/Funcion_eval(Punto_Inicial,Derdivada));

if Control_E=="Absoluto"
    Error=abs(Punto_Nuevo-Punto_Inicial);
    Error_Tabla=[0;Error];
    Punto_Inicial=Punto_Nuevo;
    Puntos=[Puntos;Punto_Inicial];
    Iteraciones=0;
    while Error>tolerancia && Iteraciones<=niter
        Punto_Nuevo=Punto_Inicial-m*(Funcion_eval(Punto_Inicial,f)/Funcion_eval(Punto_Inicial,Derdivada));
        Error=abs(Punto_Nuevo-Punto_Inicial);
        Error_Tabla=[Error_Tabla;Error];
        Punto_Inicial=Punto_Nuevo;
        Puntos=[Puntos;Punto_Inicial];
        Iteraciones=Iteraciones+1;
    end
else
    Error=abs((Punto_Nuevo-Punto_Inicial)/Punto_Nuevo);
    Error_Tabla=[0;Error];
    Punto_Inicial=Punto_Nuevo;
    Puntos=[Puntos;Punto_Inicial];
    Iteraciones=0;
    while Error>tolerancia && Iteraciones<=niter
        Punto_Nuevo=Punto_Inicial-m*(Funcion_eval(Punto_Inicial,f)/Funcion_eval(Punto_Inicial,Derdivada));
        Error=abs((Punto_Nuevo-Punto_Inicial)/Punto_Nuevo);
        Error_Tabla=[Error_Tabla;Error];
        Punto_Inicial=Punto_Nuevo;
        Puntos=[Puntos;Punto_Inicial];
        Iteraciones=Iteraciones+1;
    end
end

if Funcion_eval(Punto_Nuevo,f)==0
    Resultado=Punto_Nuevo;
    tabla=table(Puntos,Error_Tabla);
elseif Error<tolerancia
    output=sprintf('%f es una aproximación de una raiz de f(x) con una tolerancia= %f',Punto_Inicial,tolerancia);
    Resultado=Punto_Nuevo;
    tabla=table(Puntos,Error_Tabla);
else
    output=sprintf('Fracasó en %f iteraciones',niter);
end

end
