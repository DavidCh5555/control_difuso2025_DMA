%% 
%% PARAM_RHONN
% % Parámetros iniciales para la función RHONN Estructura según
% % fcn(x,u,x_hat,wk,P1k,P2k,P3k,P4k)
% 
% %% Covarianzas iniciales
P1 = eye(2)*5e2;    % para [w1 w2]
P2 = eye(2)*5e2;    % para [w3 w4]
P3 = eye(2)*5e2;    % para [w5 w6]
P4 = eye(2)*5e2;    % para [w7 w8]

%% Pesos iniciales wk = [w1 w2 w3 w4 w5 w6 w7 w8]
wk(1) = 0.04;
wk(2) = 1.1;

wk(3) = -.06;
wk(4) = 1.2;

wk(5) = 1.16;
wk(6) = -.23;

wk(7) = 1.2;
wk(8) = .04;

p1 = 5.43784e+06; p2 = 83917.5; p3 = 1295.03;

y1 = 75314.8; y2 = 9764.78; y3 = 1266.03;
% 
% %% PARAM_RHONN2-------------------------------------------------------------------------------------------------------------
% % Parámetros iniciales para la función RHONN2 Estructura según
% % fcn(x,u,x_hat,wk,P1k,P2k,P3k,P4k)
% 
%% Covarianzas iniciales
% P12 = eye(4)*5e2;    % para [w1 w2 w3 w4]
% P22 = eye(4)*5e2;    % para [w5 w6 w7 w8]
% P32 = eye(4)*5e2;    % para [w9 w10 w11 w12]
% P42 = eye(4)*5e2;    % para [w13 w14 w15 w16]
% 
% %% Pesos iniciales wk = [w1 ... w16]
% 
% wk2(1) = 0.17;
% wk2(2) = 1.03;
% wk2(3) = -.25;
% wk2(4) = -.0039;
% 
% wk2(5) = .138;
% wk2(6) = .555;
% wk2(7) = -0.17;
% wk2(8) = .43;
% 
% wk2(9) = .85;
% wk2(10) = .51;
% wk2(11) = .69;
% wk2(12) = -.63;
% 
% wk2(13) = .87;
% wk2(14) = -.32;
% wk2(15) = .88;
% wk2(16) = -.15;
% 
% p1 = 6.75327e+07;
% p2 = 509510;
% p3 = 3844.07;
% 
% y1 = 896753;
% y2 = 118243;
% y3 = 15591;

% 
% 
% % ===== Inicialización base workspace =====
% R  = 1;           % penalización escalar
% p1 = 9.84414e+10;
% p2 = 4.03853e+08;
% p3 = 1.65679e+06;
% y1 = 1.20815e+10;
% y2 = 1.58582e+09;
% y3 = 2.08156e+08;
% 
% assignin('base','R',R);
% assignin('base','p1',p1);
% assignin('base','p2',p2);
% assignin('base','p3',p3);
% assignin('base','y1',y1);
% assignin('base','y2',y2);
% assignin('base','y3',y3);
