%% model variables  (carro–péndulo 4 estados)
clear; clc;

% Símbolos
syms M m l g F
syms x1 x2 x3 x4
x = [x1; x2; x3; x4];

% Estados físicos
th    = x1;   % ángulo péndulo
dth   = x2;   % velocidad angular
x_pos = x3;   % posición carro
dx_pos= x4;   % velocidad carro

s_th = sin(th);
c_th = cos(th);

% Denominador común del modelo
den_common = m*s_th^2 + M;

% Dinámicas (mismas que en "sistema.m")
ddth   = (F*c_th + (M+m)*g*s_th - m*l*dth^2*s_th*c_th) / (den_common*l);
ddxpos = (F + m*g*s_th*c_th - m*l*dth^2*s_th) / den_common;

dx1 = dth;
dx2 = ddth;
dx3 = dx_pos;
dx4 = ddxpos;

F_vec = [dx1; dx2; dx3; dx4];

% Salida: ángulo del péndulo
G = x1;

% Jacobianos
A_var = jacobian(F_vec, x);
B_var = jacobian(F_vec, F);
C_var = jacobian(G, x);
D_var = jacobian(G, F);

% ================== Constantes numéricas ==================
M = 0.5; 
m = 0.2; 
l = 0.3; 
g = 9.81;

F  = 0;      % fuerza de equilibrio (aprox. 0)
x2 = 0;      % dθ = 0
x3 = 0;      % posición carro = 0
x4 = 0;      % velocidad carro = 0

fprintf('Linearization at -15°\n');
x1 = -15*pi/180;
A1 = eval(A_var);  B1 = eval(B_var);
C  = eval(C_var);  D  = eval(D_var);

fprintf('Linearization at -7.5°\n');
x1 = -7.5*pi/180;
A2 = eval(A_var);  B2 = eval(B_var);

fprintf('Linearization at 0°\n');
x1 = 0;
A3 = eval(A_var);  B3 = eval(B_var);

% ================== K por región (ahora 4 polos) ==================
p = [-1 -2 -3 -2000];    % AJUSTABLE: polos deseados

K_nl = place(A1, B1, p);
K_ns = place(A2, B2, p);
K_z  = place(A3, B3, p);

% Empaquetar A,B
A_all(:,:,1) = A1;  A_all(:,:,2) = A2;  A_all(:,:,3) = A3;
B_all(:,:,1) = B1;  B_all(:,:,2) = B2;  B_all(:,:,3) = B3;

% (Siguiente paso será adaptar Francis para 4x4)
 omega = 10;
 [P_nl,P_ns,P_z,G_nl,G_ns,G_z] = francis_multi_scalar(A_all,B_all,omega);
