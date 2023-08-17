clc

% PUNTO 1

disp('1)')

disp('a.')
xiyi = [8.1, 8.3, 8.7; 16.9441, 17.56492, 18.82091]';
f = Lagrange(xiyi,2) ;
f(8.4)

disp('b.')
xiyi = [0.1, 0.2, 0.3; 0.62049958, -0.25398668, 0.00660095]';
f = Lagrange(xiyi,2) ;
f(0.25)

disp('c.')
syms y
xiyi = [0, 0.5, 1, 2; 0, 25/8, 3,2]';
f = Lagrange(xiyi,3) ;
f(0.9)


%% PUNTO 4
disp('4)')
xiyi = [4,4.2,4.5,4.7,5.1,5.5,5.9,6.3,6.8,7.1;...
    102.56,113.18,130.11,142.05,167.53,195.14,224.87,256.73,299.50,326.72]' ;

disp('a.')
f = RegresionMinimosCuadrados(xiyi,1) ;
Error = sum((f(xiyi(:,1))-xiyi(:,2)).^2)

disp('b.')
f = RegresionMinimosCuadrados(xiyi,2) ;
Error = sum((f(xiyi(:,1))-xiyi(:,2)).^2)

disp('c.')
f = RegresionMinimosCuadrados(xiyi,3) ;
Error = sum((f(xiyi(:,1))-xiyi(:,2)).^2)

disp('d.')
f = RegresionMinimosCuadrados([xiyi(:,1),log(xiyi(:,2))],1) ;
y = @(x) 24.2593*exp(0.3724*x) ;
Error = sum((y(xiyi(:,1))-xiyi(:,2)).^2)
figure
plot(xiyi(:,1),xiyi(:,2),'ro',xiyi(:,1),y(xiyi(:,1)),'k')

disp('e.')
f = RegresionMinimosCuadrados([log(xiyi(:,1)),log(xiyi(:,2))],1) ;
y = @(x) exp(1.8308)*x.^(2.0195) ;
Error = sum((y(xiyi(:,1))-xiyi(:,2)).^2)
figure
plot(xiyi(:,1),xiyi(:,2),'ro',xiyi(:,1),y(xiyi(:,1)),'k')


%% PARCIAL
xiyi = [0,6,12,18,24,30,36;...
    124,134,148,152,147,133,121]' ;

N = size(xiyi,1) ;
Sx = sum(xiyi(:,1)) ;
Sy = sum(xiyi(:,2)) ;
Sxx = sum(xiyi(:,1).^2) ;
Sxxx = sum(xiyi(:,1).^3) ;
Sxxxx = sum(xiyi(:,1).^4) ;
Sxy = sum(xiyi(:,1).*xiyi(:,2)) ;
Sxxy = sum(xiyi(:,1).^2.*xiyi(:,2)) ;

A = [N,Sx,Sxx;Sx,Sxx,Sxxx;Sxx,Sxxx,Sxxxx] 
b = [Sy,Sxy,Sxxy]' ;
a = inv(A)*b
g = @(x) a(1)+a(2)*x+a(3)*x.^2 
integral(g,0,36)


%%
xiyi = [0,1,2,3;[-1,x0',1]]'
f = RegresionMinimosCuadrados(xiyi,2) ;

x1 = sum(xiyi(:,1)) ;

x2= sum(xiyi(:,1).^2) ;

x3= sum(xiyi(:,1).^3) ;

x4= sum(xiyi(:,1).^4) ;

A = [3,x1,x2;x1,x2,x3;x2,x3,x4] 

y = sum(xiyi(:,2)) ;

xy = sum(xiyi(:,1).*xiyi(:,2)) ;

x2y =sum(xiyi(:,1).^2.*xiyi(:,2)) ;

b = [y;xy;x2y]

ai = inv(A)*b

fplot(@(x) ai(3)*x^2+ai(2)*x+ai(1),[0,1])


%%
xiyi = [1,1,2;0,1,1]'
f = Lagrange(xiyi,2) ;
f = RegresionMinimosCuadrados(xiyi,1) ;
f(xiyi(:,1))


function g = RegresionMinimosCuadrados(xiyi,grado)

Terms = @(a) 0 ;
for i = 1:grado+1
     Terms = @(a) Terms(a) + a(i)*(xiyi(:,1)).^(i-1) ;
end
 
    x0 = zeros(grado+1,1)+1 ;
    ai = fsolve(@(a) Terms(a)-xiyi(:,2),x0) ;
    
syms x
P = 0 ;
for i = 1:grado+1
    P = P + ai(i)*x.^(i-1) ;
end
P
expand(P)
h = inline(P) ;

g = @(x) h(x) ;

figure('Color','White')
    fplot(g,[min(xiyi(:,1)), max(xiyi(:,1))],'k')
    hold on
    plot(xiyi(:,1),xiyi(:,2),'or')
    grid minor
    xlabel('x')
    ylabel('y = f(x)')
    

end

function g = Newton(xiyi,grado)

if grado > size(xiyi,1)-1
    error('La cantidad de datos ingresados no puede general el polinomio de grado estipulado')
end
N = grado+1 ;
a=NaN(N,N+1) ;
a(:,1)=xiyi(1,:)';
a(:,2)=xiyi(2,:)';

for i=2:N
  for j=i:N
a(j,i+1)=(a(j,i)-a(j-1,i))/(a(j,1)-a(j-i+1,1));
  end
end

 P=0 ;
for k=1:N
  dx=1 ;
  if k>1
  for i=2:k
  dx=dx*(x-xiyi(i-1,1)) ;
  end
  end
  P=P+a(k,k+1)*dx ;
end

P
f = inline(P);
g = @(x) f(x) ;
PolinomioDeNewton = simplify(P) 
    
figure('Color','White')
    fplot(g,[min(xiyi(:,1)), max(xiyi(:,1))],'k')
    hold on
    plot(xiyi(:,1),xiyi(:,2),'-or')
    grid minor
    xlabel('x')
    ylabel('y = f(x)')
    
xi=xiyi(:,1) ;
for i=2:N+1
Diferencias_divididas(:,i-1)=a(:,i) ;
end

Results = table(xi,Diferencias_divididas);
Resultados = array2table(Results)
%Gráfica
end




function g = Lagrange(xiyi,grado)

if grado > size(xiyi,1)-1
    error('La cantidad de datos ingresados no puede general el polinomio de grado estipulado')
end
N = grado+1;
syms x

P=0 ;
for k = 1:N
    
    numerador=1 ;
    denominador=1 ;
    
    for j=1:N
        if j~=k
    numerador=numerador*(x-xiyi(j,1)) ;
    denominador=denominador*(xiyi(k,1)-xiyi(j,1)) ;
    LNk = numerador/denominador ;
        end
    end
    
  P = P + xiyi(k,2)*LNk ;

end

f = inline(P);
g = @(x) f(x) ;
PolinomioDeLagrange = simplify(P) 

figure('Color','White')
    fplot(g,[min(xiyi(:,1)), max(xiyi(:,1))],'k')
    hold on
    plot(xiyi(:,1),xiyi(:,2),'or')
    grid minor
    xlabel('x')
    ylabel('y = f(x)')
end
