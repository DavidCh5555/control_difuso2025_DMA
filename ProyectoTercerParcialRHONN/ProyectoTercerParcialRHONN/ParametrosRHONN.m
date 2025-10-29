%% PARAM_RHONN_PLANTA.m
% Parámetros iniciales para la función RHONN con 4 EKF (dos pesos por estado)
% Estructura según fcn(x,u,x_hat,wk,P1k,P2k,P3k,P4k)

%% Covarianzas iniciales
P1 = eye(2)*5e2;    % para [w1 w2]
P2 = eye(2)*5e2;    % para [w3 w4]
P3 = eye(2)*5e2;    % para [w5 w6]
P4 = eye(2)*5e2;    % para [w7 w8]

%% Pesos iniciales wk = [w1 w2 w3 w4 w5 w6 w7 w8]
rng(7)   % reproducible
wk = zeros(8,1);





%% PARAM_RHONN_PLANTA.m
% Parámetros iniciales para la función RHONN con 4 EKF (cuatro pesos por estado)
% Estructura según fcn(x,u,x_hat,wk,P1k,P2k,P3k,P4k)

%% Covarianzas iniciales
P12 = eye(4)*5e2;    % para [w1 w2 w3 w4]
P22 = eye(4)*5e2;    % para [w5 w6 w7 w8]
P32 = eye(4)*5e2;    % para [w9 w10 w11 w12]
P42 = eye(4)*5e2;    % para [w13 w14 w15 w16]

%% Pesos iniciales wk = [w1 ... w16]
rng(7)   % reproducible
wk2 = zeros(16,1);

