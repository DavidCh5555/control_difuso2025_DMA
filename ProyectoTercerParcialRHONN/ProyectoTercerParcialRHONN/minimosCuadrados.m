clear; close all; clc;
disp('== IDENTIFICACIÓN LS USANDO VOLTAJES LOGGEADOS DESDE SIMULINK ==');
%% 1) Parámetros
setup_heli_2d_base;
Parametros2;
nombre_simulink = 'quanser_heli_IdentificadorRHONN';  % tu modelo
n_den = input('n_a (grado denom) para TODOS: ');
m_num = input('m_b (grado num) para TODOS: ');

%% 2) Ejecuta simulación (si ya corriste y tienes *_log en base, puedes omitir)
disp('Corriendo simulación con señales reales del modelo...');
sim(nombre_simulink);   % usa StopTime del modelo

%% 3) Leer logs desde BASE (formato esperado: matriz [t, y])
[t_vp,  Vp ] = get_log('Vp_log');     % voltaje pitch
[t_vy,  Vy ] = get_log('Vy_log');     % voltaje yaw
[t_p,   p  ] = get_log('p_log');      % ángulo pitch
[t_pdot,pd ] = get_log('pdot_log');   % vel pitch
[t_y,   y  ] = get_log('y_log');      % ángulo yaw
[t_ydot,yd ] = get_log('ydot_log');   % vel yaw

% Asegura longitudes compatibles por par
[t1, Vp,  p ] = align_pair(t_vp,  Vp,  t_p,   p );
[t2, Vp2, pd] = align_pair(t_vp,  Vp,  t_pdot,pd);
[t3, Vy,  y ] = align_pair(t_vy,  Vy,  t_y,   y );
[t4, Vy2, yd] = align_pair(t_vy,  Vy,  t_ydot,yd);

% Usa Ts del modelo si existe, si no, estima de los tiempos
if exist('Ts','var') && ~isempty(Ts)
    Ts_id = Ts;
else
    Ts_id = median(diff(t1));  % estimación
end

%% 4) Identificación por Mínimos Cuadrados
disp(' ');
disp('=========================================================================');
disp('INICIANDO ANÁLISIS DE MÍNIMOS CUADRADOS (con voltajes reales del modelo)');
disp('=========================================================================');
% Identifica los modelos (aún pueden ser inestables)
Gz_p_vp    = identificar_ls(Vp,  p,  n_den, m_num, Ts_id, '1. Pitch Angle (p / Vp)');
Gz_pdot_vp = identificar_ls(Vp2, pd, n_den, m_num, Ts_id, '2. Pitch Velocity (p\_dot / Vp)');
Gz_y_vy    = identificar_ls(Vy,  y,  n_den, m_num, Ts_id, '3. Yaw Angle (y / Vy)');
Gz_ydot_vy = identificar_ls(Vy2, yd, n_den, m_num, Ts_id, '4. Yaw Velocity (y\_dot / Vy)');

%% 4b) ESTABILIZACIÓN DE MODELOS (NUEVA SECCIÓN) %% <-- CAMBIO
% Revisa cada modelo y fuerza la estabilidad si es necesario.
disp(' ');
disp('=========================================================================');
disp('VERIFICANDO ESTABILIDAD Y ESTABILIZANDO MODELOS');
disp('=========================================================================');
Gz_p_vp    = estabilizar_modelo(Gz_p_vp,    '1. Gz(p / Vp)');
Gz_pdot_vp = estabilizar_modelo(Gz_pdot_vp, '2. Gz(p_dot / Vp)');
Gz_y_vy    = estabilizar_modelo(Gz_y_vy,    '3. Gz(y / Vy)');
Gz_ydot_vy = estabilizar_modelo(Gz_ydot_vy, '4. Gz(y_dot / Vy)');

%% 5) Resumen
disp(' ');
disp('=========================================================================');
disp('RESUMEN DE FUNCIONES DE TRANSFERENCIA (FINALES Y ESTABLES)');
disp('=========================================================================');
disp('Modelo 1: Gz(p / Vp)');     disp(Gz_p_vp);
disp('Modelo 2: Gz(p_dot / Vp)'); disp(Gz_pdot_vp);
disp('Modelo 3: Gz(y / Vy)');     disp(Gz_y_vy);
disp('Modelo 4: Gz(y_dot / Vy)'); disp(Gz_ydot_vy);
disp('== FIN ==');

%% 6) Estados iniciales (para "Initial states" de cada bloque)
% Esto ahora usa los modelos ESTABILIZADOS de la sección 4b
[x0_p_vp]    = x0_from_data(Gz_p_vp,    Vp,  p);
[x0_pdot_vp] = x0_from_data(Gz_pdot_vp, Vp2, pd);
[x0_y_vy]    = x0_from_data(Gz_y_vy,    Vy,  y);
[x0_ydot_vy] = x0_from_data(Gz_ydot_vy, Vy2, yd);
disp(' '); disp('================ ESTADOS INICIALES (ESTABLES, copia tal cual) ================');
print_x0('p/Vp    ', x0_p_vp);
print_x0('p_dot/Vp', x0_pdot_vp);
print_x0('y/Vy    ', x0_y_vy);
print_x0('y_dot/Vy', x0_ydot_vy);

function print_x0(name,x0)
    fprintf('%s : order=%d  x0 = %s\n', name, numel(x0), mat2str(x0(:).', 8));
end

%% ==================== Funciones auxiliares NUEVAS ====================

function x0 = x0_from_data(Gz, u, y)
    % Estados iniciales consistentes con filtro en z^-1
    [b,a] = tfdata(Gz,'v');          % coef. tal cual
    na = numel(a)-1; nb = numel(b)-1;
    ord = max(na,nb);
    if numel(y) <= ord || numel(u) <= ord
        x0 = zeros(ord,1); return
    end
    y_hist = flipud(y(end-ord:end-1));   % y(k-1..k-ord)
    u_hist = flipud(u(end-ord:end-1));   % u(k-1..k-ord)
    yic = zeros(1,na);  uic = zeros(1,nb);
    yic(1:min(na,ord)) = y_hist(1:min(na,ord));
    uic(1:min(nb,ord)) = u_hist(1:min(nb,ord));
    x0 = filtic(b,a,yic,uic);            % vector columna
end

%% ==================== FUNCIÓN DE ESTABILIZACIÓN (NUEVA) ==================== %% <-- CAMBIO
function Gz_stable = estabilizar_modelo(Gz_inestable, nombre_modelo)
    % Revisa polos. Si alguno es >= 1.0, lo mueve a 0.9999
    polos = pole(Gz_inestable);
    idx_inestables = (abs(polos) >= 1.0);
    
    if ~any(idx_inestables)
        % El modelo ya es estable
        fprintf('  %s: Modelo ya es estable. (Polos max: %f)\n', nombre_modelo, max(abs(polos)));
        Gz_stable = Gz_inestable;
        return;
    end
    
    % El modelo es inestable, hay que corregir
    fprintf('! %s: ¡Modelo INESTABLE detectado! Polos: \n', nombre_modelo);
    disp(polos(idx_inestables).');
    
    % Obtener componentes del modelo
    ceros = zero(Gz_inestable);
    ganancia = Gz_inestable.K;
    Ts = Gz_inestable.Ts;
    
    % Corregir polos: moverlos de >=1.0 a 0.9999
    polos_corregidos = polos;
    polos_corregidos(idx_inestables) = 0.9999 * (polos_corregidos(idx_inestables) ./ abs(polos_corregidos(idx_inestables)));
    
    fprintf('  Polos corregidos a: \n');
    disp(polos_corregidos(idx_inestables).');
    
    % Reconstruir el modelo
    Gz_stable = zpk(ceros, polos_corregidos, ganancia, Ts);
    % Convertir de nuevo al formato TF con variable z^-1
    Gz_stable = tf(Gz_stable, 'variable', 'z^-1');
end

%% ==================== Funciones auxiliares ====================

function [t, y] = get_log(varname)
    % Lee variable del base workspace y la normaliza a [t,y]
    if ~evalin('base', sprintf('exist(''%s'',''var'')', varname))
        error('No existe %s en el workspace. Revisa tus MATLAB Function logger.', varname);
    end
    L = evalin('base', varname);
    if isnumeric(L) && size(L,2) >= 2
        t = L(:,1);
        y = L(:,2);
    elseif istimetable(L)
        t = seconds(L.Time - L.Time(1));
        y = L.Variables(:,1);
    elseif isstruct(L) && isfield(L,'time') && isfield(L,'signals')
        t = L.time(:);
        y = L.signals.values(:);
    else
        error('Formato no soportado para %s. Esperado: matriz [t,y], timetable, o structure-with-time.', varname);
    end
    m = ~(isnan(t)|isnan(y));
    t = t(m); y = y(m);
end

function [t, u2, y2] = align_pair(tu, u, ty, y)
    % Alinea por índice al mínimo largo. Asume mismo Ts o casi.
    N = min(numel(u), numel(y));
    u = u(:); y = y(:);
    tu = tu(:); ty = ty(:);
    u2 = u(1:N);
    y2 = y(1:N);
    t  = (0:N-1) * median([median(diff(tu)), median(diff(ty))], 'omitnan');
end

function Gz_hat = identificar_ls(u, y, n_den, m_num, Ts, fig_title)
    % --- Sección original de LS ---
    u = u(:); y = y(:);
    N = length(y);
    if length(u) ~= N, error(['U e Y no coinciden: ' fig_title]); end
    L = n_den + m_num;
    Phi = zeros(N, L);
    Y   = y;
    for k = 1:N
        y_pas = zeros(1,n_den);
        for i = 1:n_den
            if k-i >= 1, y_pas(i) = -y(k-i); end
        end
        u_pas = zeros(1,m_num);
        for j = 1:m_num
            if k-j >= 1, u_pas(j) =  u(k-j); end
        end
        Phi(k,:) = [y_pas, u_pas];
    end
    Theta = (Phi' * Phi) \ (Phi' * Y);
    a_hat = Theta(1:n_den);
    b_hat = Theta(n_den+1:end);
    Gz_hat = tf(b_hat', [1 a_hat'], Ts, 'variable','z^-1');
    
    % --- Gráficas (MODIFICADAS) --- %% <-- CAMBIO
    
    % 1. Predicción a 1-paso (la que se ve bien)
    y_hat = Phi * Theta;
    e = Y - y_hat;
    t = (0:N-1)*Ts;
    
    % 2. Simulación 'free-run' (la que muestra la inestabilidad)
    %    Se simula con x0=0 solo para diagnóstico visual de estabilidad
    y_sim = lsim(Gz_hat, u, t); 
    
    figure('Name', fig_title);
    
    % Gráfica comparativa
    subplot(2,1,1);
    plot(t,Y,'b', t,y_hat,'r--', t,y_sim,'g-.', 'LineWidth',1.3);
    title(['Comparación - ' fig_title]); grid on;
    legend('y (Real)', 'y^ (Predicción 1-paso)', 'y_{sim} (Simulación x0=0)', 'Location', 'best');
    xlabel('t [s]'); ylabel('ampl');
    % Ajusta límites Y si la simulación explota, para poder ver lo demás
    ylim_sim = [min(Y), max(Y)];
    if any(abs(y_sim) > 1e4) % Si explota
        ylim_sim = [min(Y)*1.5, max(Y)*1.5];
    end
    ylim(ylim_sim);
    
    % Gráfica de error (de la predicción)
    subplot(2,1,2);
    plot(t,e,'k','LineWidth',1.0); grid on; xlabel('t [s]'); ylabel('e');
    title('Error de predicción (e = y - y^)');
end