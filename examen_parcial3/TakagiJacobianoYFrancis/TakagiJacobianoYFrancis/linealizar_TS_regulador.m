%% linearizar_TS_regulador.m
% Recalcula K_nl, K_ns, K_z, K_ps, K_pl para tu controlador TS puro.
% Usa la misma planta (pendulo) y los mismos puntos del programador B: [-15 -7.5 0 7.5 15] grados.

% --- Constantes de la planta (ajusta si tu 'sistema.m' cambia) ---
l = 1; m = 1; g = 9.81; b = 0.5;
I = m*l^2;

% --- Linealizaciones en grados: [-15 -7.5 0 7.5 15] ---
angles_deg = [-15 -7.5 0 7.5 15];
A_all = zeros(2,2,5);
B_all = zeros(2,1,5);
for i = 1:5
    th = deg2rad(angles_deg(i));
    A_all(:,:,i) = [0 1; (g/l)*cos(th)  -(b*l)/I];
    B_all(:,:,i) = [0; 1/I];
end

% --- Polos deseados (par complejo conjugado): ajusta zeta y wn si quieres ---
zeta = 0.7; wn = 6;                    % amortiguamiento y rapidez
wd = wn*sqrt(max(0,1 - zeta^2));
p   = [-zeta*wn + 1j*wd, -zeta*wn - 1j*wd];   % 2 polos

% --- Ganancias locales por punto ---
K_nl = place(A_all(:,:,1), B_all(:,:,1), p);   % -15°
K_ns = place(A_all(:,:,2), B_all(:,:,2), p);   % -7.5°
K_z  = place(A_all(:,:,3), B_all(:,:,3), p);   % 0°
K_ps = place(A_all(:,:,4), B_all(:,:,4), p);   % 7.5°
K_pl = place(A_all(:,:,5), B_all(:,:,5), p);   % 15°

% --- Publica al Workspace como filas 1x2 (lo que espera tu Control.m) ---
K_nl = K_nl(:).';  K_ns = K_ns(:).';  K_z = K_z(:).';
K_ps = K_ps(:).';  K_pl = K_pl(:).';
assignin('base','K_nl',K_nl);
assignin('base','K_ns',K_ns);
assignin('base','K_z', K_z);
assignin('base','K_ps',K_ps);
assignin('base','K_pl',K_pl);

fprintf('TS listo: K_nl=%s, K_ns=%s, K_z=%s, K_ps=%s, K_pl=%s\n', ...
    mat2str(K_nl,4), mat2str(K_ns,4), mat2str(K_z,4), mat2str(K_ps,4), mat2str(K_pl,4));
