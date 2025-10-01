function [P, G] = EcuacionesDeFrancisDMA(A, B, C, D, S, R)
% FRANCIS_SOLVER_MEJORADO
%
% Resuelve las Ecuaciones de Francis (Eigenstructure Assignment) para encontrar
% las matrices de transformación P y ganancia G.
% El sistema resuelto es:
% 1. A*P + B*G = P*S
% 2. C*P + D*G = R
%
% Entradas:
%   A, B, C, D: Matrices del modelo de estado (x_dot = Ax + Bu, y = Cx + Du).
%   S: Matriz de autovalores/modos deseados (r x r).
%   R: Vector/matriz de condiciones de salida deseadas (p x r).
%
% Salidas:
%   P: Matriz de transformación (n x r).
%   G: Matriz de ganancia (m x r).
%

    % 1. Obtener dimensiones del sistema
    [n, ~] = size(A);    % n = Número de estados
    [~, m] = size(B);    % m = Número de entradas
    [p, ~] = size(C);    % p = Número de salidas/restricciones
    [r, ~] = size(S);    % r = Número de autovalores deseados (tamaño del subsistema)

    % 2. Construcción de las sub-matrices del sistema lineal M*z = rhs
    %    Se usa la identidad kron(X, Y) * vec(Z) = vec(Y Z X.')
    %    El vector de incógnitas es z = [vec(P); vec(G.')]

    % --- Bloque superior (derivado de A*P + B*G = P*S) ---
    % Ecuación de Lyapunov extendida: A*P - P*S = -B*G
    % Forma matricial linealizada: (kron(I, A) - kron(S.', I)) * vec(P) + kron(I, B) * vec(G) = 0
    
    % Matriz para vec(P)
    M1_P = kron(eye(r), A) - kron(S.', eye(n));   % Dimensiones: (n*r) x (n*r)
    
    % Matriz para vec(G)
    M1_G = kron(eye(r), B);                      % Dimensiones: (n*r) x (m*r)
    
    % --- Bloque inferior (derivado de C*P + D*G = R) ---
    % Forma matricial linealizada: kron(I, C) * vec(P) + kron(I, D) * vec(G) = vec(R)
    
    % Matriz para vec(P)
    M2_P = kron(eye(r), C);                      % Dimensiones: (p*r) x (n*r)
    
    % Matriz para vec(G): Manejo del caso D=0 (si no se proporciona)
    if isempty(D) || isequal(D, 0)
        M2_G = zeros(p * r, m * r);
    else
        M2_G = kron(eye(r), D);                  % Dimensiones: (p*r) x (m*r)
    end
    
    % --- Vector Lado Derecho (RHS) ---
    % El término R debe ir al lado derecho y con signo negativo en el bloque 2.
    rhs = [zeros(n * r, 1); -R(:)]; % R(:) vectoriza la matriz R (p*r x 1)

    % 3. Ensamblaje de la matriz de coeficientes global M
    M = [M1_P, M1_G;
         M2_P, M2_G];
    
    % 4. Resolver el sistema lineal
    % z = M \ rhs;  (z = [vec(P); vec(G.)])
    z = M \ rhs;

    % 5. Reestructurar las soluciones a sus formas matriciales P y G
    % P (n x r): La matriz P se vectoriza en el primer bloque de z.
    P_vector = z(1 : n * r);
    P = reshape(P_vector, [n, r]);

    % G (m x r): La matriz G se vectoriza en el segundo bloque de z.
    G_vector = z(n * r + 1 : end);
    % Nota: La vectorización usa vec(G.') en la formulación, por lo que reshape(G_vector, [m, r]) daría G.'.
    % Hacemos reshape a [r, m] y transponemos para obtener G (m x r).
    G_transpuesto = reshape(G_vector, [r, m]);
    G = G_transpuesto.';

end