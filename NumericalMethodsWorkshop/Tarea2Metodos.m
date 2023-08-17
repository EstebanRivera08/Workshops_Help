clear
clc

%%

f = @(x) 9.*x./(1-3.*x) ;
x0 = 0 ;
grado = 4;
p = 0.1 ;
Taylor(f,x0,grado,p)


%% Punto 1
f = @(x) x.^3 ;
x0 = 0 ;
grado = 2;
p = 0.5 ;
Taylor(f,x0,grado,p)

%Punto 2
f = @(x) (x-1).*log(x) ;
x0 = 1 ;
grado = 3;
p = 0.5 ;
Taylor(f,x0,grado,p)

%Punto 3
f = @(x) x.*exp(x.^2) ;
x0 = 0 ;
grado = 4;
p = 0.4 ;
Taylor(f,x0,grado,p)

%%

%Punto 4
disp('Punto 4')

disp('a)')
f = @(x) x-2.^(-x) ;

a = 0 ; b = 1 ;
Boltzano(f,a,b,10^-5)

%%
disp('b)')
f = @(x) exp(x)-x.^2+3*x-2 ;

a = 0 ; b = 1 ;
Boltzano(f,a,b,10^-5)

%%

disp('c)')
f = @(x) 2*x.*cos(2*x)-(x+1).^2 ;

a = -3 ; b = -2 ;
Boltzano(f,a,b,10^-5)

a = -1 ; b = 0 ;
Boltzano(f,a,b,10^-5)

%%

disp('d)')
f = @(x) -2*x.^2+3*x-1 ;

% a = 0.2 ; b = 0.3 ;
% Boltzano(f,a,b,10^-5)

% a = 1.2 ; b = 1.3 ;
% Boltzano(f,a,b,10^-5)
hold on
fplot(f,[0,1.5],'k')
plot([0.5,1],[0,0],'or')
yline(0,'b')
ylabel('y = f(x)')
xlabel('x')
grid minor

%%
disp('Punto 5')

disp('a)')

f = @(x) x.^3 -2*x.^2-5 ;

a = 1 ; b = 4 ;
NR(f,(a+b)/2,10^-4)


%%

disp('b)')
f = @(x) x.^3 + 3*x.^2-1 ;

a = -3 ; b = -2 ;
NR(f,(a+b)/2,10^-4)

%%
disp('c)')
f = @(x) x-cos(x) ;

a = 0 ; b = pi/2 ;
NR(f,(a+b)/2,10^-4)

%%
disp('d)')
f = @(x) x-0.8*sin(x) ;

a = 0 ; b = pi/2 ;
NR(f,(a+b)/2,10^-4)



function g = Taylor(f,x0,grado,p)

T = f(x0) ;

syms x
for n=1:grado
    df = inline(diff(f(x),n)) ;
    T = T + df(x0)*(x-x0).^n/factorial(n) ;
end

disp('Polinomio de Taylor:')
T 
dfn_1 = inline(diff(f(x),grado+1)) ;

disp('Función de error del polinomio de Taylor:')
R = dfn_1(p)*(x-x0)^(grado+1)/factorial(grado+1) 
Taylor = inline(T) ;
Error = inline(R) ;

disp('Resultados de la aproximación en el punto p:')
table(f(p),Taylor(p),abs(f(p)-Taylor(p)),Error(p),'VariableNames',{'f(x)','P(x)','P(x)-f(x)','R(x)'})

g = Taylor ;
end



function g = NR(f,x0,Error)
    IterMax = 100 ;
    syms x
    
    df=inline(diff(f(x))) ;

    if abs(df(x0))<1e-10
      error('La función no es diferenciable en el punto ingresado')         
    end        

    
%Condiciones iniciales
    k=0 ; 
%Iteraciones
   Resultados(k+1,1)=k ; Resultados(k+1,2)=x0 ;
   Resultados(k+1,3)=f(x0) ; Resultados(k+1,4)=df(x0) ;
   xc = 1000 ;
   
while ((abs(xc-x0))>Error && k<IterMax)
   
   %Newton Raphson
   xc = x0 ;
   x=x0-f(x0)/df(x0) ; 
   x0=x ;
   k=k+1 ;    %#iteración
          %Respuesta cada iteración 
   Resultados(k+1,1)=k ; Resultados(k+1,2)=x0 ;
   Resultados(k+1,3)=f(x0) ; Resultados(k+1,4)=df(x0) ;
   
end  

%Código gráfica
    figure
    y=[x0-10:0.1:x0+10] ; %limites de la gráfica
    hold on
    yline(0,'b')
    plot(y,f(y),'k',Resultados(k,2),Resultados(k,3), 'ro')
    ylim([-10 10])
    ylabel('y = f(x)')
    xlabel('x')
    grid minor

%Codigo para el formato de la tabla
format short
NR=array2table(Resultados,...
            'VariableNames',{'i','x(i)','f(x(i))','df(xi))'})    
fprintf('La cantidad de iteraciones fueron : %f \n', k)          
fprintf('La raiz de la función en el intervalo dado es x = : %f \n', x0)

g = x0 ;
end



function g = Boltzano(f,a,b,Error)

    IterMax = 100 ;
    if f(a)*f(b)>0
       error('No existe una raiz en el intervalo estipulado, por favor ingrese otros limites') 
    end     
                 
    k=1 ;
    c=(a+b)/2 ; %punto inicial
    
   Resultados(1,1)=0 ; Resultados(1,2)=a ;  Resultados(k+1,3)=c ;  
   Resultados(1,4)=b ; Resultados(1,5)=f(a) ;
   Resultados(1,7)=f(c) ; Resultados(1,6)=f(b) ;      
   Err = Error + 1 ;
  
while Err>=Error && k<IterMax  
    
    % Respuestas cada iteración
   Resultados(k,1)=k ; Resultados(k+1,2)=a ;  Resultados(k+1,3)=c ;  
   Resultados(k,4)=b ; Resultados(k+1,5)=f(a) ;
   Resultados(k,7)=f(c) ; Resultados(k+1,6)=f(b) ;  
    
       if (f(a)*f(c)>0)        
         a=c ;
        end 
        
        if (f(a)*f(c)<0)        
         b=c ;        
        end
        
   c=(a+b)/2 ; % Redefinir C
   k=k+1 ;    %#iteración
   Err=abs(b-a) ;  %Tamaño intervalo
                   
end  

%Código gráfica
    figure
    x=c-10:0.1:c+10 ; %limites de la gráfica
    hold on
    yline(0,'b')
    plot(x,f(x),'k',Resultados(k-1,3),Resultados(k-1,7), 'or')
    axis tight
    ylim([-10 10])
    ylabel('y = f(x)')
    xlabel('x')
    grid minor
    

%Codigo para el formato de la tabla
format short
BB=array2table(Resultados(1:k-1,:),...
            'VariableNames',{'i','a','c','b','f(a)','f(c)','f(b)'})    
fprintf('La cantidad de iteraciones fue : %f \n', k-1)        
fprintf('La raiz de la función en el intervalo dado es x = : %f \n', c)        
       
g = c ;
end    






