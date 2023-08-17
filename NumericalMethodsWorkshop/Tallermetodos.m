format long
clear all
clc



%%
disp('PUNTO 1')

f = @(x) 1/sqrt(2*pi)*exp(-1/2*x.^2) ;

disp('a)')

a = -1 ; b = 1 ; N = 50 ;
fprintf("\n La probabilidad de que esté entre [-sigma,sigma], es: \n")
pm = punto_medio_comp(f,a,b,N) ;
simp = simpson_comp(f,a,b,N)  ;
trap = trapecio_comp(f,a,b,N)  ;
fprintf(' Punto medio: %0.6f \n Simpson Compuesto: %0.6f \n Trapecio compuesto: %0.6f \n\n',pm,simp,trap)
Errorpm = abs(integral(f,a,b)-pm)
Errorsimp = abs(integral(f,a,b)-simp)
Errortrap = abs(integral(f,a,b)-trap)
disp('b)')

a = -2 ; b = 2 ; N = 50 ;
fprintf("\n La probabilidad de que esté entre [-2sigma,2sigma], es: \n")
pm = punto_medio_comp(f,a,b,N) ;
simp = simpson_comp(f,a,b,N)  ;
trap = trapecio_comp(f,a,b,N)  ;
fprintf(' Punto medio: %0.6f \n Simpson Compuesto: %0.6f \n Trapecio compuesto: %0.6f \n\n',pm,simp,trap)
Errorpm = abs(integral(f,a,b)-pm)
Errorsimp = abs(integral(f,a,b)-simp)
Errortrap = abs(integral(f,a,b)-trap)
disp('c)')

a = -3 ; b = 3 ; N = 50 ;
fprintf("\n La probabilidad de que esté entre [-3sigma,3sigma], es: \n")
pm = punto_medio_comp(f,a,b,N) ;
simp = simpson_comp(f,a,b,N)  ;
trap = trapecio_comp(f,a,b,N)  ;
fprintf(' Punto medio: %0.6f \n Simpson Compuesto: %0.6f \n Trapecio compuesto: %0.6f \n\n',pm,simp,trap)
Errorpm = abs(integral(f,a,b)-pm)
Errorsimp = abs(integral(f,a,b)-simp)
Errortrap = abs(integral(f,a,b)-trap)

disp('PUNTO 2')
Tiempo = [0,6,12,18,24,30,36,42,48,54,60,66,72,78,84] ;
Velocidad = [124,134,148,156,147,133,121,109,99,85,78,89,104,116,123] ;
N = size(Tiempo,2)-1 ;
h = (Tiempo(end)-Tiempo(1))/N ;
Impar = 0 ; Par = 0 ;
    for i=2:N
        if ((-1)^i)==1
           Impar=Impar+Velocidad(i) ;
        else
           Par=Par+Velocidad(i) ;
        end
    end
Distancia = h/3*(Velocidad(1)+4*Impar+2*Par+Velocidad(end)) ;
fprintf('Empleando el método de simpson encontramos que la \ndistancia recorrida será: %.2f pies \n',Distancia)

disp('PUNTO 3')
f1 = @(x) sqrt(x).*log(x) ;
f2 = @(x) exp(-x).*sin(x)./x ;
f3 = @(x) exp(-x.^4) ;

%Nodos
m = [2 3 5 7 10 15] ;
%Particiones
N = m-1 ;
for i=1:size(m,2)
    Result(i,1) = cuadratura_gaussiana(f1,0,1,N(i)) ;
    Result(i,2) = cuadratura_gaussiana(f2,0,inf,N(i)) ;
    Result(i,3) = cuadratura_gaussiana(f3,-inf,inf,N(i)) ;
end

figure('Color','White')
plot(m,Result)
xlabel('No. Nodos')
ylabel('Valor integral')
legend('\surdx ln(x)')
grid minor

function g = punto_medio_comp(f,a,b,N)
    h = (b-a)/N ;
    x = linspace(a,b,N+1) ;
    x_prom = (x(2:end)+x(1:end-1))/2 ;
    g = h*sum(f(x_prom)) ;
end

function g = simpson_comp(f,a,b,N)
    h = (b-a)/N ;
    x = linspace(a,b,N+1) ;
    Impar = 0 ; Par = 0 ;
    for i = 2:N
        if ((-1)^i)==1
           Impar=Impar+f(x(i)) ;
        else
           Par=Par+f(x(i)) ;
        end
    end
    g = h/3*(f(a)+4*Impar+2*Par+f(b)) ;
end

function g = trapecio_comp(f,a,b,N)
    h = (b-a)/N ;
    x = linspace(a,b,N+1) ;
    g = h/2*(f(a)+f(b)+2*sum(f(x(2:N)))) ;
end
    
function g = CuadGauss2p(f,a,b)
    x0 = -1/sqrt(3) ;
    x1 = -x0 ;
    xa = (b+a)/2 + (b-a)/2*x0 ;
    xb = (b+a)/2 + (b-a)/2*x1 ;
    g = ((b-a)/2)*(f(xa) + f(xb)) ;
end

function g = cuadratura_gaussiana(f,a,b,N)
    x = linspace(a,b,N+1)' ;
    area = 0 ;
        %Para las gráficas
    x0 = -1/sqrt(3) ;
    x1 = -x0 ;
    for i=1:N
        deltaA =  CuadGauss2p(f,x(i),x(i+1)) ;
        area = area + deltaA ;   
    end
    
    g = area ;
end