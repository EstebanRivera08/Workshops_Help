
clc
clear
%Planta ubicada en santander

%INFORMACIÓN

%Datos del proceso
TGin = 23; % °C, Temp. bulbo seco = Temp. promedio Santander
Twin = TGin-5; % °C
TLin = 45 ; % °C
TLout = 25 ; %°C
L = 15 ; %kg/s

%Rapidez mínima
P = 1 ;  %atm
Lmin = 2.7 ; %kg/m2s
Gmin = 1.8 ; %kg/m2s

%Dureza en las corrientes
xC = 2000 ;
xM = 500 ;

%CÁLCULOS PRELIMINARES
Rango = TLin-TLout
Aproximacion = TLout-Twin
Efectividad = Rango/(Aproximacion+Rango)

% COMPARACIÓN
rhoAire = 1.107 ;%kg/m3
CpAire = 1007 ; %J/kgK
kAire = 0.0279 ; %W/mK
muAire = 1.872e-5 ; %kg/ms
DAB = 2.58e-5 ; %m2/s
Sc = muAire/rhoAire/DAB 
Pr = muAire*CpAire/kAire
Le = Sc/Pr
RelLewis = Le^0.567

%PROPIEDADES DEL AIRE Y AGUA

%Masa molecular
MA = 18.02 ; %kg/kmol
MAS = 28.67 ; %kg/kmol

%Capacidades calorificas
CpAire = 1.005 ; %kJ/kgAS°C
CpVap = 1.88 ; %kJ/kgA°C
CpAL = 4.187 ; %kJ/kgA°C

%Capacidad calorífica humeda
Cph = @(Y) CpAire + Y*CpVap ; %kJ/kgAS
T0 = 0 ;

%Presión de saturación
Psat = @(T) 10.^(8.07131 - 1730.630./(233.426+T))/760 ; %atm

%Calor de vaporización
C1 = 5.2053e7 ; C2 = 0.3199 ; C3 =-0.212 ; C4 = 0.25795; 
Tc = 647.14 ; %K
lambda =  @(T) (10^-3*C1.*(1-(T+273.15)/Tc).^(C2+C3.*(T+273.15)/Tc+C4.*((T+273.15)/Tc).^2))/MA ; %kJ/kg

%CALCULOS PSICROMETRICOS
%Humedad
Ys = @(T,P) MA/MAS*(Psat(T)./(P-Psat(T))) ; %kgA/kgAS
YG = @(TG,Tw,P) Ys(Tw,P) - Cph(Ys(Tw,P)).*(TG-Tw)/lambda(Tw) ; %kgA/kgAS

%Entalpia
HG = @(TG,Tw,P) (CpAire+YG(TG,Tw,P)*CpVap)*(TG-T0)+YG(TG,Tw,P)*lambda(T0) ; %kJ/kgAS°C
HGs = @(T,P) (CpAire+Ys(T,P)*CpVap).*(T-T0)+Ys(T,P)*lambda(T0) ; %kJ/kgAS°C

%CURVAS
HGin = HG(TGin,Twin,P)
YGin = YG(TGin,Twin,P)

%Entalpía de equilibrio y su derivada
Hs = @(T) HGs(T,P) ;
syms x
dHs = eval(['@(x)' char(diff(Hs(x)))])  ;

%Curva de operación
C_Opemin = @(T,H) (HGin-H)/(TLout-TLin)*(T-TLout)+HGin ;

Obj = @(T,H) [((HGin-H)/(TLout-TLin)-dHs(T)); C_Opemin(T,H)-Hs(T)]  ;
X = fsolve(@(X) Obj(X(1),X(2)),[37;200]) ;
Ttoque = X(1) ;

%Entalpia de salida mínima
HGout_min = X(2) 

Gsmin = L*CpAL*(TLin-TLout)/(HGout_min-HGin)

Gs = 1.5*Gsmin 
L

%MODELOS DEL EMPAQUE

% Área
A = min([L/Lmin,Gs/Gmin])
% A = pi*D^2/4
D = sqrt(A/pi*4) ;

G_A = Gs/(1+YGin)/A
L_A = L/A 

Kya = 1.98*(L_A)^0.45*(G_A)^0.75  % kg/sm2
DP =  27.8*(L_A)^0.35*(G_A)^0.55 % atm

%Curva de operación real
H_Ope =@(T) L/Gs*CpAL*(T-TLout)+HGin ;
HGout = H_Ope(TLin)
YGout = Ys(Ttoque,P)
TGout = Ttoque

NtOG1 = integral(@(T) 1./(Hs(T)-H_Ope(T)),TLout,TLin) 

Tint = TLout:0.01:TLin ;
for i = 1:size(Tint,2)
    df(i) = dHs(Tint(i)) ;
end

NtOG2 = trapz(Hs(Tint), 1./(Hs(Tint)-H_Ope(Tint)))
%NtOG3 = trapz(H_Ope(Tint),1./(Hs(Tint)-H_Ope(Tint)))

HtOG1 = L*CpAL/Kya/A
HtOG2 = L/A/Kya
%HtOG3 = Gs/A/Kya

V1 = NtOG1*L*CpAL/Kya
V2 = NtOG2*L/Kya
Z1 = NtOG1*HtOG1
Z2 = NtOG2*HtOG2
%Z3 = NtOG3*HtOG3

[z,TL] = ode45(@(z,TL) Kya*A*(Hs(TL)-H_Ope(TL))/dHs(TL)/L,[0 Z2],TLout) ;

% Cálculo Make-up

E = Gs*(YGout-YGin)
W = 0.2/100*Gs
B = E*xM/(xC-xM)-W
M = B+E+W

% Q = F*R*T/P = (kg/s)*(kmol/kg)*(kPa*m3/kmol/K)*(K)/(kPa)= m3/s
QAire = Gs/28.97*(8.3145)*(TGin+273.15)/(P*101.325)*3600

hp = DP/9.8  % mm H2O

QAgua = L/1000  %m3/s

DT = 2.54*1.5 ; %cm
AT = pi/4*(DT/100)^2 % m2
VAgua = QAgua/AT

muAgua = 1e-3 ; %Pas
rhoAgua = 1000; %kg/m3

Re = rhoAgua*DT/100*VAgua/muAgua

f = fsolve(@(f) 1/sqrt(f)+2*log(2.51/(Re*sqrt(f))),0.01)

Lt = 4 ; %m

hf = 4*f*Lt/(DT/100)*VAgua^2/2/9.8


%% GRAFICAS
T = 20:50 ;
figure('Color','White')
hold on
% 'color',[0.9290 0.6940 0.1250]
plot(T,HGs(T,P),'k')
plot([TLin,TLout],[HGout_min,HGin],'bo-')
plot([TLin,TLout],[H_Ope(TLin),H_Ope(TLout)],'ro-')
plot(Ttoque,Hs(Ttoque),'sk')
xline(TLin,'--')
xlabel('Temperatura del liquido, C','interpreter','latex')
ylabel('Entalpia Aire-Vapor agua, kJ/kgAS','interpreter','latex')
ylim([50,250])
grid minor
legend('H_{Eq}','H_{Ope. Min.}','H_{Ope. Real}','location','northwest')

figure('Color','White')
tiledlayout(1,2,'tilespacing','compact','padding','compact')

nexttile
hold on
% 'color',[0.9290 0.6940 0.1250]
plot(T,1./(Hs(T)-H_Ope(T)),'k')
fill([TLout,Tint,TLin],[0,1./(Hs(Tint)-H_Ope(Tint)),0],'b')
xlabel('Temperatura del liquido, C','interpreter','latex')
ylabel('$\frac{1}{H_{Eq}-H_{G}^{Ope.}}, \frac{kgAS}{kJ}$ ','interpreter','latex','fontsize',16)
grid minor

nexttile
hold on
% 'color',[0.9290 0.6940 0.1250]
plot(Hs(T),1./(Hs(T)-H_Ope(T)),'k')
fill([Hs(TLout),Hs(Tint),Hs(TLin)],[0,1./(Hs(Tint)-H_Ope(Tint)),0],'r')
xlabel('$ H_{Eq}, \frac{kJ}{kgAS} $','interpreter','latex')
ylabel('$\frac{1}{H_{Eq}-H_{G}^{Ope.}}, \frac{kgAS}{kJ}$ ','interpreter','latex','fontsize',16)
xlim([50,250])
grid minor

figure('color','white')
plot([0 Z2],[TGin, TGout],'b',z,TL,'r',[0 Z2],[TLout, TLin],'r--')
xlabel('$ Z_{torre}, m $','interpreter','latex')
ylabel('$Temperatura, C$','interpreter','latex')
xline(Z2,'--')
legend('T_G','T_L')
grid minor

figure('color','white')
plot([0 Z1],[YGin, YGout],'k--',[0 Z2],[YGin, YGout],'k')
xlabel('$ Z_{torre}, m $','interpreter','latex')
ylabel('$Humedad absoluta, \frac{kgH_2O}{kgAS}$','interpreter','latex')
legend('Z_1','Z_2')
grid minor





