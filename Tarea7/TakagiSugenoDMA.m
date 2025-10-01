%% Furuta: linealizaciones TS, ganancias K_i y (P_i,G_i) de Francis
clear; clc;
%% Parámetros (paper)
m2 = 0.50;     l2 = 0.75;   l1 = 0.12;   J  = 0.003;
Km = 0.104;    Kg = 0.055;  Rm = 1.9;    g  = 9.81;
P1 = (Kg*Km)/Rm;
P2 = (Kg^2*Km^2)/Rm;

%% Símbolos
syms th al dth dal v real
x = [th; al; dth; dal];

% Matrices D,C,G (del modelo dinámico)
d11 = m2*l1^2 + J;
d12 = m2*l1*l2*cos(al);
d21 = d12;
d22 = m2*l2^2;
c11 = 0;
c12 = -m2*l1*l2*sin(al)*dal;
c21 =  m2*l1*l2*sin(al)*dth;
c22 = 0;
g1 = 0;
g2 =  m2*l2*g*sin(al);

D = [d11 d12; d21 d22];
C_mat = [c11 c12; c21 c22]; % Renombrado a C_mat para no confundir con C de la salida
G_vec = [g1; g2];           % Renombrado a G_vec

qdot = [dth; dal];
tau  = P1*v - P2*dth;
rhs  = [tau; 0] - C_mat*qdot - G_vec;
qdd  = D \ rhs;

f = [dth;
     dal;
     qdd(1);
     qdd(2)];

g_in = diff(f, v);          % Vector de entrada (∂f/∂v)
f    = subs(f, v, 0);       % f(x, 0)
A_var = jacobian(f, x);     % Matriz A (∂f/∂x)
B_var = g_in;               % Matriz B

%% Salida y exosistema de Francis
C = [1 0 0 0];              % Restricción de salida: y = theta
D_francis = 0;              % Matriz D de la restricción de salida (escalar)
W = 1;                  % Frecuencia del exosistema (autovalores: +/- i*omega)
S = [0 W; -W 0];    % Matriz de autovalores deseados para el subsistema de Francis (r=2)
R = [1 0];                  % Condición de salida deseada (R, p x r): y_desired = R * z_francis

%% Puntos TS (premisa = alpha en grados)
alphas_deg = [-15 -7.5 0 7.5 15];
alphas = alphas_deg*pi/180;
A_all = zeros(4,4,5);
B_all = zeros(4,1,5);
for k = 1:5
    % th=0, dth=0, dal=0 en los puntos de equilibrio
    A_all(:,:,k) = double(subs(A_var, {th,al,dth,dal}, {0, alphas(k), 0, 0}));
    B_all(:,:,k) = double(subs(B_var, {al}, {alphas(k)}));
end

%% Ganancias LQR/polos por región (ejemplo con pole placement)
% Este código permanece igual, simplemente calcula las ganancias K para los 5 modelos
p = [-.1 -3 -8 -2];   % Polos deseados para el subespacio controlable
K = zeros(1,4,5);
I4 = eye(4);
for k = 1:5
    A = A_all(:,:,k);  B = B_all(:,:,k);
    [Abar,Bbar,~,T,kidx] = ctrbf(A,B,I4);   % Transformación a forma canónica de controlabilidad
    nc = sum(kidx);                          % dimensión del subespacio controlable
    A11 = Abar(1:nc,1:nc);
    B1  = Bbar(1:nc,:);
    
    % Coloca polos solo en la parte controlable (usa tantos polos como nc)
    Kc = place(A11,B1,p(1:nc));              % 1×nc
    
    % Mapea al sistema original: K = [Kc  0] * T^{-1}
    K_full = [Kc, zeros(1, size(A,1)-nc)] / T;
    K(1,1:4,k) = K_full;
end

K_nl = squeeze(K(:,:,1)); % Nivel Bajo Negativo
K_ns = squeeze(K(:,:,2)); % Nivel Medio Negativo
K_z  = squeeze(K(:,:,3)); % Zero (Vertical)
K_ps = squeeze(K(:,:,4)); % Nivel Medio Positivo
K_pl = squeeze(K(:,:,5)); % Nivel Alto Positivo


%% Francis por región
% Aquí se llama a la función EcuacionesDeFrancisDMA (tu solucionador general)
% para cada punto de operación (modelo TS).

% Modelo 1: alpha = -15 deg
[P_nl,G_nl] = EcuacionesDeFrancisDMA(A_all(:,:,1), B_all(:,:,1), C, D_francis, S, R);

% Modelo 2: alpha = -7.5 deg
[P_ns,G_ns] = EcuacionesDeFrancisDMA(A_all(:,:,2), B_all(:,:,2), C, D_francis, S, R);

% Modelo 3: alpha = 0 deg
[P_z ,G_z ] = EcuacionesDeFrancisDMA(A_all(:,:,3), B_all(:,:,3), C, D_francis, S, R);

% Modelo 4: alpha = 7.5 deg
[P_ps,G_ps] = EcuacionesDeFrancisDMA(A_all(:,:,4), B_all(:,:,4), C, D_francis, S, R);

% Modelo 5: alpha = 15 deg
[P_pl,G_pl] = EcuacionesDeFrancisDMA(A_all(:,:,5), B_all(:,:,5), C, D_francis, S, R);

% Opcional: Mostrar algunos resultados para verificar
disp('--- Matrices P y G de Francis para el punto de equilibrio CERO (alpha=0) ---');
disp('P_z = '); disp(P_z);
disp('G_z = '); disp(G_z);