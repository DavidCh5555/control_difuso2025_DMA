function [P_nl,P_ns,P_z,G_nl,G_ns,G_z] = francis_multi_scalar(A_all,B_all,omega)
% FRANCIS_MULTI_SCALAR
%   Diseña reguladores tipo Francis para cada región.
%   A_all(:,:,i): matrices A (4x4)
%   B_all(:,:,i): matrices B (4x1)
%   omega: frecuencia del exosistema (dw=[0 w;-w 0]w)
%
%   Devuelve:
%   P_* : matrices P (4x2)
%   G_* : matrices G (1x2)

    [P_nl,G_nl] = solve_one(A_all(:,:,1),B_all(:,:,1),omega);
    [P_ns,G_ns] = solve_one(A_all(:,:,2),B_all(:,:,2),omega);
    [P_z ,G_z ] = solve_one(A_all(:,:,3),B_all(:,:,3),omega);
end

% =========================================================
function [P,G] = solve_one(A,B,omega)
% Resuelve las ecuaciones de Francis:
%   A P + B G = P S
%   C P       = I_1x2
%
% con:
%   S = [0  ω; -ω 0]
%   C = [1 0 0 0]  (seguimos regulando el ángulo del péndulo x1)
%
% A: nxn, B: nx1  -> P: nx2, G: 1x2

    n = size(A,1);
    S = [0 omega; -omega 0];
    C = [1 zeros(1,n-1)];   % salida: primer estado (ángulo)

    I2 = eye(2);
    In = eye(n);

    % Ecuación 1: A P + B G = P S
    % (I2⊗A - S.'⊗In) vec(P) + (I2⊗B) vec(G) = 0
    M1   = [ kron(I2,A) - kron(S.',In),  kron(I2,B) ];
    rhs1 = zeros(2*n,1);

    % Ecuación 2: C P = [1 0]
    % (I2⊗C) vec(P) = [1; 0]
    M2   = [ kron(I2,C), zeros(2,2) ];
    rhs2 = [1; 0];

    % Sistema lineal total
    M   = [M1; M2];
    rhs = [rhs1; rhs2];

    sol  = M \ rhs;
    vecP = sol(1:n*2);
    vecG = sol(n*2+1:end);

    P = reshape(vecP, [n, 2]);   % nx2
    G = vecG.';                  % 1x2
end
