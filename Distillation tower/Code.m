clc
clear

load Datos

%% Diagramas de equilibrio

P = 1.01325 ; %bar

PsatT = @(T) 10.^(4.05004 - 1327.62./(217.625+T)) ; %T(°C) P(bar)
PsatCH = @(T) 10.^(3.93002 - 1182.774./(220.618+T)) ; %T(°C) P(bar)

TebCH = 80.75 ; %°C
TebT = 110.6; %°C

T = linspace(TebCH, TebT,50) ;
x1 = (P-PsatT(T))./(PsatCH(T)-PsatT(T)) ;
y1 = PsatCH(T)/P.*x1 ;

figure('Color','white')
tiledlayout(1,2,'padding','compact','tilespacing','compact')
nexttile
plot(x1,y1,'b',[0 1],[0 1],'k')
xlabel('x_{1}')
ylabel('y_{1}')
axis([0 1 0 1])
grid minor
legend('Equilibrio','y = x','location','southeast')

nexttile
plot(x1,T,'b',y1,T,'r')
xlabel('x_{1}, y _{1}')
ylabel('T , °C')
xlim([0 1])
grid minor
legend('Temp. de burbuja','Temp. de rocío')

%% Cálculo de Rmin

% Información de la columna de destilación

%Alimentación
F = 150 ; %mol/h
zF_CH = 0.45 ;
zF_T = 0.55 ;

%Destilado
xD_CH = .96 ;
xD_T = 1-xD_CH ;

%Fondos
xB_T = .96 ;
xB_CH = 1- xB_T ;

% Creamos curva de equilibrio
y_Eq = @(x) interp1(x1,y1,x,'spline') ;

%Valor de q sub enfriado
q = 1.5 ;

%Recta del balance en el plato de alimentación
y_int = @(x) q/(q-1)*x - 1/(q-1)*zF_CH ;

%Intersección RBPA y la curva de equilibrio
x_BPA = fsolve(@(x) y_Eq(x)-y_int(x),zF_CH)

% Pendiente curva de operación R/(R+1)
m = (y_int(x_BPA)-xD_CH)/(x_BPA-xD_CH)

%Cálculo de Rmin
Rm = m/(1-m)

figure('Color','white')
hold on
plot(x1,y1,'r',[0 1],[0 1],'k')
plot([x_BPA,zF_CH],y_int([x_BPA,zF_CH]),'color',[0.9290 0.6940 0.1250])
plot([xB_CH,x_BPA,xD_CH],[xB_CH,y_int(x_BPA),xD_CH],'b-s')
xlabel('x_{1}')
ylabel('y_{1}')
axis([0 1 0 1])
grid minor
legend('Equilibrio','y = x','Balance en el plato de alimentación','Recta de operación mínima','location','southeast')


%% Cálculo de la curva de operación y de las etapas

alfa = 1.4 ;

%Cálculo de relación de relujo de operación
R = alfa*Rm

%Recta de operación en la zona de enriquecimiento
y_Ope1 = @(x) R/(1+R)*x+xD_CH/(1+R) ;

%Intersección y_Ope1 y y_Eq
x_int = fsolve(@(x) y_Ope1(x)-y_int(x),zF_CH)

%Recta de operación en la zona de despojamiento
y_Ope2 = @(x) interp1([x_int,xB_CH],[y_int(x_int),xB_CH],x,'linear') ;

%Creamos la función de xEq = f(yEq)
x_Eq = @(x) interp1(y1,x1,x,'spline') ;

%Composiciones en las etapas de equilibrio
x_etapas(1) = xD_CH ;
y_etapas(1) = xD_CH ;
k = 1;
i = 1;
while x_etapas(i) > xB_CH
    x_etapas(i+1) = x_Eq(y_etapas(i)) ;
    y_etapas(i+1) = y_Ope1(x_etapas(i+1));
    if x_etapas(i+1) > x_int
        x_etapas(i+1) = x_Eq(y_etapas(i)) ;
        y_etapas(i+1) = y_Ope1(x_etapas(i+1));
    else 
        e2(k) = i ;
        x_etapas(i+1) = x_Eq(y_etapas(i));
        y_etapas(i+1) = y_Ope2(x_etapas(i+1)) ;     
        k = k+1;
    end
    i = i+1 ;
end

%Reorganización de las etapas y No de etapas
NoEtapas = i - 1
EtapaFeed = e2(1) 
x_Etapas = sort([x_etapas,x_etapas],'descend') ;
y_Etapas = sort([y_etapas,y_etapas],'descend') ;

figure('Color','white')
hold on
plot(x1,y1,'r',[0 1],[0 1],'k')
plot([x_BPA,zF_CH],y_int([x_BPA,zF_CH]),'color',[0.9290 0.6940 0.1250])
plot([xB_CH,x_int,xD_CH],[xB_CH,y_int(x_int),xD_CH],'b-s')
plot(x_Etapas(1:end-1),y_Etapas(2:end),'k')
xlabel('x_{1}')
ylabel('y_{1}')
axis([0 1 0 1])
grid minor
legend('Equilibrio','y = x','Balance en el plato de alimentación',...
    'Recta de operación','Etapas','location','southeast')

%% Etapas mínimas de equilibrio

%Composiciones en las etapas de equilibrio
x_min(1) = xD_CH ;
y_min(1) = xD_CH ;

i = 1;
k = 1 ;
while x_min(i) > xB_CH && i < 50
    x_min(i+1) = x_Eq(y_min(i)) ;
    y_min(i+1) = x_min(i+1);
    if x_min(i+1) > x_int
        x_min(i+1) = x_Eq(y_min(i)) ;
        y_min(i+1) = x_min(i+1);
    else 
        e2(k) = i ;
        x_min(i+1) = x_Eq(y_min(i));
        y_min(i+1) = x_min(i+1) ;  
        k = k+1 ;
    end
    i = i+1 ;
end

%Reorganización de las etapas y No de etapas
NoEtapas_min = i - 1
EtapaFeed_min = e2(1) 
x_min = sort([x_min,x_min],'descend') ;
y_min = sort([y_min,y_min],'descend') ;

figure('Color','white')
hold on
plot(x1,y1,'r',[0 1],[0 1],'k')
plot([x_BPA,zF_CH],y_int([x_BPA,zF_CH]),'color',[0.9290 0.6940 0.1250])
plot([xB_CH,x_int,xD_CH],[xB_CH,y_int(x_int),xD_CH],'b-s')
plot(x_min(2:end),y_min(1:end-1),'k')
xlabel('x_{1}')
ylabel('y_{1}')
axis([0 1 0 1])
grid minor
legend('Equilibrio','y = x','Balance en el plato de alimentación',...
    'Recta de operación','Etapas mínimas','location','southeast')
%% Cargas de condensador y Rehervidor

%Temperatura crítica
Tc_CH = 553.80 ;%K
Tc_T = 591.75 ; %K

%Constantes de la función de calor de vaporización
C1 = [4.4902e7,0.39881,0,0] ;
C2 = [4.9507e7,0.37742,0,0] ;  

%Capacidades caloríficas promedio del gas y del líquido
CpL_CH = (1.4836+2.0323)/2*1e5*1e-3 ; %kj/kmolK 
CpL_T = (2.3774+1.3507)/2*1e5*1e-3 ; %kj/kmolK
CpV_CH = 1.125*1e5*1e-3 ; %kJ/kmolK
CpV_T = 1.195*1e5*1e-3 ; %kJ/kmolK

%Funciones de entalpía
Tref = 25 ; %°C
Hv = @(y1,T) y1*(Hvap(C1,Tc_CH,Tref)+CpV_CH*(T-Tref))+...
    (1-y1)*(Hvap(C2,Tc_T,Tref)+CpV_T*(T-Tref)) ; %J/molK

Hl = @(x1,T) x1*(CpL_CH*(T-Tref))+(1-x1)*(CpL_T*(T-Tref)) ; %J/molK

%% Simulación

% Temperatura de burbuja y de rocío
%Burbuja
Tx = @(x) interp1(x1,T,x,'spline') ;

%Rocío
Ty = @(x) interp1(y1,T,x,'spline') ;

%Cálculo temperatura de entrada
TF = fsolve(@(Tf) q - (Hl(zF_CH,Tf)-Hv(zF_CH,Ty(zF_CH)))...
    ./(Hl(zF_CH,Tx(zF_CH))-Hv(zF_CH,Ty(zF_CH))),50)

%Light and heavy key recov
% F = D + B
% F*zF = D*xD+B*xB = D*xD + (F-D)*xB
% -> D = F*(zF-xB)/(zF-xD)

% Flujos molares
D = F*(zF_CH-xB_CH)/(xD_CH-xB_CH)
B = F-D

L = R*D
V = L + D
L_ = L + q*F
V_ = V + (q-1)*F

Qreb = V_*Hv(y_etapas(end-1),Ty(y_etapas(end-1))) +...
    B*Hl(xB_CH,Tx(xB_CH))-L_*Hl(x_etapas(end),Tx(x_etapas(end)))

Qcond = V*Hv(y_etapas(1),Ty(y_etapas(1))) - ...
    D*Hl(xD_CH,Tx(xD_CH))-L*Hl(x_etapas(1),Tx(x_etapas(1)))

Distillatetofeedratio = D/F
lightkeyrecov = D*xD_CH/(F*zF_CH)
heavykeyrecov = D*(1-xD_CH)/(F*(1-zF_CH))

%Relación de reflujo
R
Rm

%Número de etapas
NoEtapas
EtapaFeed


%% Comparación métodos termodinámicos

figure('Color','white')
hold on
p1 = plot(1-x1,T+273.15,'r',1-y1,T+273.15,'r') ;
p2 = plot(VLE1(:,3),VLE1(:,1),'ko',VLE1(:,5),VLE1(:,1),'ko');
p3 = plot(VLE3(:,3),VLE3(:,1),'k^',VLE3(:,5),VLE3(:,1),'k^');
p4 = plot(VLE5(:,3),VLE5(:,1),'ks',VLE5(:,5),VLE5(:,1),'ks');
p5 = plot(VLE7(:,3),VLE7(:,1),'kd',VLE7(:,5),VLE7(:,1),'kd');
p6 = plot(SRK(:,1),SRK(:,2),'b',SRK(:,3),SRK(:,2),'b');
p7 = plot(PSRK(:,1),PSRK(:,2),'k',PSRK(:,3),PSRK(:,2),'k');
xlabel('x_{Tol}, y _{Tol}')
ylabel('T , °C')
grid minor
xlim([0 1])
legend([p2(1),p3(1),p4(1),p5(1),p1(1),p6(1),p7(1)],'Exp. Data [1]','Exp. Data [2]',...
    'Exp. Data [3]','Exp. Data [4]','Solución ideal','SRK','PSRK','location','northwest')

function f=Hvap(Ci,Tc,T)
    T = T+273.15 ;
    f = 10^-3*Ci(1)*(1-(T/Tc)).^(Ci(2)+Ci(3)*(T/Tc)+Ci(4)*(T/Tc).^2) ; %kJ/kmol
end



