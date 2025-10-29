%% minimosCuadradosMIMO.m
clear; close all; clc;
disp('== IDENTIFICACIÓN LS MISO (MIMO) USANDO VOLTAJES LOGGEADOS ==');

%% 1) Parámetros
setup_heli_2d_base;
Parametros2;
nombre_simulink = 'quanser_heli_IdentificadorRHONN';  % tu modelo
n_den = input('n_a (grado denom común) para TODOS: ');
m_num = input('m_b (grado num entradas) para TODOS: ');

%% 2) Ejecuta simulación (si ya corriste y tienes *_log en base, puedes omitir)
disp('Corriendo simulación con señales reales del modelo...');
sim(nombre_simulink);   % usa StopTime del modelo

%% 3) Leer logs desde BASE (formato esperado: matriz [t, y])
[t_vp,  Vp ] = get_log('Vp_log');     % Input 1: voltaje pitch
[t_vy,  Vy ] = get_log('Vy_log');     % Input 2: voltaje yaw
[t_p,   p  ] = get_log('p_log');      % Output 1: ángulo pitch
[t_pdot,pd ] = get_log('pdot_log');   % Output 2: vel pitch
[t_y,   y  ] = get_log('y_log');      % Output 3: ángulo yaw
[t_ydot,yd ] = get_log('ydot_log');   % Output 4: vel yaw

% Asegura tipos y Ts
Vp = double(Vp(:)); Vy = double(Vy(:));
p  = double(p(:));  pd = double(pd(:));
y  = double(y(:));  yd = double(yd(:));

if exist('Ts','var') && ~isempty(Ts), Ts_id = double(Ts);
else, Ts_id = median(diff(double(t_vp))); end

%% 3b) ALINEACIÓN GLOBAL
N_min = min([numel(Vp), numel(Vy), numel(p), numel(pd), numel(y), numel(yd)]);
disp(['Recortando todas las señales a N = ' num2str(N_min) ' muestras.']);
Vp = Vp(1:N_min); Vy = Vy(1:N_min);
p  = p(1:N_min);  pd = pd(1:N_min);
y  = y(1:N_min);  yd = yd(1:N_min);
t  = (0:N_min-1)' * Ts_id;

%% 4) Identificación por Mínimos Cuadrados MISO
disp(' ');
disp('=========================================================================');
disp('INICIANDO ANÁLISIS DE MÍNIMOS CUADRADOS MISO (2 Entradas, 1 Salida)');
disp('=========================================================================');

[Gz_p_vp,  Gz_p_vy ] = identificar_ls_miso(Vp, Vy, p,  n_den, m_num, Ts_id, '1. p | (Vp,Vy)');
[Gz_pd_vp, Gz_pd_vy] = identificar_ls_miso(Vp, Vy, pd, n_den, m_num, Ts_id, '2. p_dot | (Vp,Vy)');
[Gz_y_vp,  Gz_y_vy ] = identificar_ls_miso(Vp, Vy, y,  n_den, m_num, Ts_id, '3. y | (Vp,Vy)');
[Gz_yd_vp, Gz_yd_vy] = identificar_ls_miso(Vp, Vy, yd, n_den, m_num, Ts_id, '4. y_dot | (Vp,Vy)');

%% 4b) ESTABILIZACIÓN
disp(' ');
disp('=========================================================================');
disp('VERIFICANDO ESTABILIDAD Y ESTABILIZANDO MODELOS');
disp('=========================================================================');
Gz_p_vp   = estabilizar_modelo(Gz_p_vp,   'Gz(p/Vp)');
Gz_p_vy   = estabilizar_modelo(Gz_p_vy,   'Gz(p/Vy)');
Gz_pd_vp  = estabilizar_modelo(Gz_pd_vp,  'Gz(pd/Vp)');
Gz_pd_vy  = estabilizar_modelo(Gz_pd_vy,  'Gz(pd/Vy)');
Gz_y_vp   = estabilizar_modelo(Gz_y_vp,   'Gz(y/Vp)');
Gz_y_vy   = estabilizar_modelo(Gz_y_vy,   'Gz(y/Vy)');
Gz_yd_vp  = estabilizar_modelo(Gz_yd_vp,  'Gz(yd/Vp)');
Gz_yd_vy  = estabilizar_modelo(Gz_yd_vy,  'Gz(yd/Vy)');

%% 5) Resumen
disp(' ');
disp('=========================================================================');
disp('RESUMEN DE FUNCIONES DE TRANSFERENCIA MISO (FINALES Y ESTABLES)');
disp('=========================================================================');
disp('Modelo 1: p(z)     = Gz(p/Vp)*Vp(z)   + Gz(p/Vy)*Vy(z)');
disp('Modelo 2: p_dot(z) = Gz(pd/Vp)*Vp(z)  + Gz(pd/Vy)*Vy(z)');
disp('Modelo 3: y(z)     = Gz(y/Vp)*Vp(z)   + Gz(y/Vy)*Vy(z)');
disp('Modelo 4: y_dot(z) = Gz(yd/Vp)*Vp(z)  + Gz(yd/Vy)*Vy(z)');

%% 6) Estados iniciales para MATLAB Function (x0 por rama)
[x01_p,  x02_p ] = x0_from_data_miso(Gz_p_vp,  Gz_p_vy,  Vp, Vy, p);
[x01_pd, x02_pd] = x0_from_data_miso(Gz_pd_vp, Gz_pd_vy, Vp, Vy, pd);
[x01_y,  x02_y ] = x0_from_data_miso(Gz_y_vp,  Gz_y_vy,  Vp, Vy, y);
[x01_yd, x02_yd] = x0_from_data_miso(Gz_yd_vp, Gz_yd_vy, Vp, Vy, yd);

disp(' ');
disp('================ ESTADOS INICIALES MISO (copia tal cual) ================');
print_x0('Salida p     | x01 (Vp→p)',  x01_p);
print_x0('Salida p     | x02 (Vy→p)',  x02_p);
print_x0('Salida p_dot | x01 (Vp→pd)', x01_pd);
print_x0('Salida p_dot | x02 (Vy→pd)', x02_pd);
print_x0('Salida y     | x01 (Vp→y)',  x01_y);
print_x0('Salida y     | x02 (Vy→y)',  x02_y);
print_x0('Salida y_dot | x01 (Vp→yd)', x01_yd);
print_x0('Salida y_dot | x02 (Vy→yd)', x02_yd);

%% 7) BLOQUES PARA PEGAR EN MATLAB Function (coeficientes + x0)
disp(' ');
disp('=========================================================================');
disp('=============== INICIO DE DATOS PARA MATLAB FUNCTION ====================');
disp('=========================================================================');
fprintf('Ts = %g;\n\n', Ts_id);
print_coeffs_for_fcn('PITCH (p)',           Gz_p_vp,  Gz_p_vy,  x01_p,  x02_p,  n_den, m_num);
print_coeffs_for_fcn('PITCH_VEL (pd)',      Gz_pd_vp, Gz_pd_vy, x01_pd, x02_pd, n_den, m_num);
print_coeffs_for_fcn('YAW (y)',             Gz_y_vp,  Gz_y_vy,  x01_y,  x02_y,  n_den, m_num);
print_coeffs_for_fcn('YAW_VEL (yd)',        Gz_yd_vp, Gz_yd_vy, x01_yd, x02_yd, n_den, m_num);
disp('=========================================================================');
disp('================ FIN DE DATOS PARA MATLAB FUNCTION ======================');
disp('=========================================================================');

%% ==================== FUNCIONES AUXILIARES ====================

function [Gz_u1_hat, Gz_u2_hat] = identificar_ls_miso(u1, u2, y, n_den, m_num, Ts, fig_title)
    u1 = u1(:); u2 = u2(:); y = y(:);
    N = length(y);
    if length(u1) ~= N || length(u2) ~= N
        error(['Longitudes no coinciden en: ' fig_title]); 
    end
    L = n_den + m_num + m_num;
    Phi = zeros(N, L);
    Y   = y;

    for k = 1:N
        y_pas = zeros(1, n_den);
        for i = 1:n_den
            if k-i >= 1, y_pas(i) = -y(k-i); end
        end
        u1_pas = zeros(1, m_num);
        for j = 1:m_num
            if k-j >= 1, u1_pas(j) =  u1(k-j); end
        end
        u2_pas = zeros(1, m_num);
        for j = 1:m_num
            if k-j >= 1, u2_pas(j) =  u2(k-j); end
        end
        Phi(k,:) = [y_pas, u1_pas, u2_pas];
    end

    % Regularización ligera para robustez
    lambda = 5e-1;
    Theta = (Phi' * Phi + lambda*eye(L)) \ (Phi' * Y);

    idx = 1;
    a_hat  = Theta(idx : idx+n_den-1);         idx = idx + n_den;
    b1_hat = Theta(idx : idx+m_num-1);         idx = idx + m_num;
    b2_hat = Theta(idx : end);

    den_comun = [1 a_hat'];
    Gz_u1_hat = tf(b1_hat', den_comun, Ts, 'variable','z^-1');
    Gz_u2_hat = tf(b2_hat', den_comun, Ts, 'variable','z^-1');

    % Predicción 1-paso y simulación libre
    y_hat = Phi * Theta;
    t = (0:N-1)*Ts;
    y_sim = lsim(Gz_u1_hat, u1, t) + lsim(Gz_u2_hat, u2, t);

    figure('Name', fig_title);
    subplot(2,1,1);
    plot(t,Y,'b', t,y_hat,'r--', t,y_sim,'g-.', 'LineWidth',1.2);
    title(['Comparación MISO - ' fig_title]); grid on;
    legend('y (Real)', 'y^ (1-paso)', 'y_{sim} (x0=0)', 'Location', 'best');
    xlabel('t [s]'); ylabel('ampl');

    subplot(2,1,2);
    plot(t, Y - y_hat, 'k', 'LineWidth',1.0); grid on; xlabel('t [s]'); ylabel('e');
    title('Error de predicción (e = y - y^)');
end

function Gz_stable = estabilizar_modelo(Gz_in, nombre_modelo)
    polos = pole(Gz_in);
    idx_bad = (abs(polos) >= 1.0);
    if ~any(idx_bad)
        fprintf('  %s: Modelo ya es estable. (Polos max: %g)\n', nombre_modelo, max(abs(polos)));
        Gz_stable = Gz_in; return;
    end
    fprintf('! %s: ¡Modelo INESTABLE detectado! Polos: \n', nombre_modelo);
    disp(polos(idx_bad).');
    ceros = zero(Gz_in); Ts = Gz_in.Ts; zpk_model = zpk(Gz_in);
    K = zpk_model.K;
    polos(idx_bad) = 0.9999 * (polos(idx_bad)./abs(polos(idx_bad)));
    fprintf('  Polos corregidos a: \n'); disp(polos(idx_bad).');
    Gz_stable = tf(zpk(ceros, polos, K, Ts)); Gz_stable.Variable = 'z^-1';
end

function [x01, x02] = x0_from_data_miso(Gz_u1, Gz_u2, u1, u2, y)
    % x0 para dos ramas SISO en paralelo (z^-1), mismo denominador a.
    u1 = double(u1(:)); u2 = double(u2(:)); y = double(y(:));
    [b1,a] = tfdata(Gz_u1,'v');  b1 = double(b1);  a = double(a);
    [b2,~] = tfdata(Gz_u2,'v');  b2 = double(b2);
    na = numel(a)-1; ord = na;
    if numel(y) <= ord || numel(u1) <= ord || numel(u2) <= ord
        x01 = zeros(ord,1); x02 = zeros(ord,1); return
    end
    % Históricos (k-1..k-ord)
    y_hist  = flipud(y(end-ord:end-1));
    u1_hist = flipud(u1(end-ord:end-1));
    u2_hist = flipud(u2(end-ord:end-1));

    % Aproxima y2_hist pasando u2 por G2 con IC=0 (recursión ARX)
    y2_hist = zeros(ord,1);
    nb2 = numel(b2)-1;
    for k=1:ord
        acc = 0;
        for j=1:min(k,nb2)
            acc = acc + b2(j+1)*u2_hist(k-j+1); % b2(j+1) multipl. u2(k-j)
        end
        for i=1:min(k,na)
            acc = acc - a(i+1)*y2_hist(k-i+1);
        end
        y2_hist(k) = acc;
    end
    y1_hist = y_hist - y2_hist;

    % Estados tipo DF1 por rama con FILTIC
    nb1 = numel(b1)-1;
    yic1 = zeros(1,na); uic1 = zeros(1,nb1);
    yic2 = zeros(1,na); uic2 = zeros(1,nb2);

    yic1(1:min(na,ord)) = y1_hist(1:min(na,ord)).';
    uic1(1:min(nb1,ord)) = u1_hist(1:min(nb1,ord)).';
    yic2(1:min(na,ord)) = y2_hist(1:min(na,ord)).';
    uic2(1:min(nb2,ord)) = u2_hist(1:min(nb2,ord)).';

    x01 = filtic(b1, a, yic1, uic1);
    x02 = filtic(b2, a, yic2, uic2);
end

function [t, y] = get_log(varname)
    if ~evalin('base', sprintf('exist(''%s'',''var'')', varname))
        error('No existe %s en el workspace. Revisa tus MATLAB Function logger.', varname);
    end
    L = evalin('base', varname);
    if isnumeric(L) && size(L,2) >= 2
        t = double(L(:,1)); y = double(L(:,2));
    elseif istimetable(L)
        t = seconds(L.Time - L.Time(1)); y = double(L.Variables(:,1));
    elseif isstruct(L) && isfield(L,'time') && isfield(L,'signals')
        t = double(L.time(:)); y = double(L.signals.values(:));
    else
        error('Formato no soportado para %s. Esperado: matriz [t,y], timetable, o structure-with-time.', varname);
    end
    m = ~(isnan(t)|isnan(y)); t = t(m); y = y(m);
end

function print_x0(name,x0)
    fprintf('%s : order=%d  x0 = %s\n', name, numel(x0), mat2str(x0(:).', 10));
end

function print_coeffs_for_fcn(model_name, Gz_u1, Gz_u2, x01, x02, n_den, m_num)
    [b1,a] = tfdata(Gz_u1,'v'); [b2,~] = tfdata(Gz_u2,'v');
    a_full  = [a,  zeros(1, n_den + 1 - numel(a))];
    b1_full = [b1, zeros(1, m_num + 1 - numel(b1))];
    b2_full = [b2, zeros(1, m_num + 1 - numel(b2))];
    x01_full = [x01(:).', zeros(1, numel(a_full)-1 - numel(x01))];
    x02_full = [x02(:).', zeros(1, numel(a_full)-1 - numel(x02))];

    fprintf('--- %s ---\n', model_name);
    fprintf('a   = %s;\n',   mat2str(a_full, 16));
    fprintf('b1  = %s;  %% u1→y\n', mat2str(b1_full, 16));
    fprintf('b2  = %s;  %% u2→y\n', mat2str(b2_full, 16));
    fprintf('x01 = %s;  %% x0 rama u1\n', mat2str(x01_full, 16));
    fprintf('x02 = %s;  %% x0 rama u2\n\n');
end
