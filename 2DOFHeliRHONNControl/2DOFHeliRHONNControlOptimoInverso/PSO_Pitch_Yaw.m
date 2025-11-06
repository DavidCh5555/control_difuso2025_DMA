%% =============== CONFIG ===============
mdl        = 'HeliRHONNControl';   % tu .slx (sin extensión)
max_iter   = 10000;    % iteraciones
pop_size   = 50;       % particulas
c1         = 0.20;     % factor cognitivo (0..1)
c2         = 0.30;     % factor social    (0..1)
partial_log = true;    % imprime progreso por iteración

% Warm-start manual
use_init   = true;
init_pitch = [9.84414e+10 4.03853e+08 1.65679e+06];
init_yaw   = [1.20815e+10 1.58582e+09 2.08156e+08];

% Para depurar solo pitch (ignora yaw en el coste)
ignore_yaw = false;

% usar RMSE ya calculado en Simulink / workspace ===
use_rmse_from_workspace = true;          % << activa para NO calcular MSE aquí
rmse_pitch_name = 'RMSE_pitch';          % nombre de la señal/variable de RMSE de pitch
rmse_yaw_name   = 'RMSE_yaw';            % nombre de la señal/variable de RMSE de yaw
% Si tus bloques se llaman distinto, cambia estos nombres arriba.

% Semilla (fija para reproducibilidad, comenta para aleatorio)
rng(1);

% Asegura que el modelo esté cargado y actualizado
load_system(mdl);
set_param(mdl,'ReturnWorkspaceOutputs','on','SimulationMode','normal');
set_param(mdl,'FastRestart','off');
set_param(mdl,'SimulationCommand','update');

% Estructura partícula (3 variables por eje)
empty_particle.position = [];
empty_particle.velocity = [];
empty_particle.cost     = inf;
empty_particle.best.position = [];
empty_particle.best.cost     = inf;

% Enjambres separados
pitch = repmat(empty_particle, pop_size, 1);
yaw   = repmat(empty_particle, pop_size, 1);

% Global best por eje
gb_pitch = empty_particle.best; gb_pitch.cost = inf;
gb_yaw   = empty_particle.best; gb_yaw.cost   = inf;

% ===== Inicialización =====
for i = 1:pop_size
    if use_init
        if i==1
            pitch(i).position = init_pitch(:);
            yaw(i).position   = init_yaw(:);
        else
            pitch(i).position = init_pitch(:).* (1 + 0.05*(randn(3,1)));
            yaw(i).position   = init_yaw(:)  .* (1 + 0.05*(randn(3,1)));
        end
    else
        pitch(i).position = [1+10*rand; 10*randn; 1+10*rand];
        yaw(i).position   = [1+10*rand; 10*randn; 1+10*rand];
    end
    pitch(i).velocity = zeros(3,1);
    yaw(i).velocity   = zeros(3,1);

    % Proyecta a SPD (ajusta hasta que eig>0)
    pitch(i).position = project_to_spd_vec(pitch(i).position);
    yaw(i).position   = project_to_spd_vec(yaw(i).position);

    % Evalúa coste inicial (ahora puede leer RMSE directo)
    [Jp, Jy] = eval_cost(mdl, pitch(i).position, yaw(i).position, ignore_yaw, ...
                         use_rmse_from_workspace, rmse_pitch_name, rmse_yaw_name);
    pitch(i).cost = Jp; yaw(i).cost = Jy;
    pitch(i).best = struct('position', pitch(i).position, 'cost', pitch(i).cost);
    yaw(i).best   = struct('position', yaw(i).position,   'cost', yaw(i).cost);

    if pitch(i).best.cost < gb_pitch.cost, gb_pitch = pitch(i).best; end
    if yaw(i).best.cost   < gb_yaw.cost,   gb_yaw   = yaw(i).best;   end
end

best_costs_pitch = nan(max_iter,1);
best_costs_yaw   = nan(max_iter,1);

% ===== Bucle PSO =====
for it = 1:max_iter
    for j = 1:pop_size
        % Velocidad y posición (Pitch)
        r1 = rand(3,1); r2 = rand(3,1);
        pitch(j).velocity = pitch(j).velocity ...
            + c1.*r1.*(pitch(j).best.position - pitch(j).position) ...
            + c2.*r2.*(gb_pitch.position     - pitch(j).position);
        pitch(j).position = project_to_spd_vec(pitch(j).position + pitch(j).velocity);

        % Velocidad y posición (Yaw)
        r1 = rand(3,1); r2 = rand(3,1);
        yaw(j).velocity = yaw(j).velocity ...
            + c1.*r1.*(yaw(j).best.position - yaw(j).position) ...
            + c2.*r2.*(gb_yaw.position     - yaw(j).position);
        yaw(j).position = project_to_spd_vec(yaw(j).position + yaw(j).velocity);

        % Evalúa
        [Jp, Jy] = eval_cost(mdl, pitch(j).position, yaw(j).position, ignore_yaw, ...
                             use_rmse_from_workspace, rmse_pitch_name, rmse_yaw_name);
        pitch(j).cost = Jp; yaw(j).cost = Jy;

        % Mejores locales
        if pitch(j).cost < pitch(j).best.cost
            pitch(j).best.position = pitch(j).position;
            pitch(j).best.cost     = pitch(j).cost;
        end
        if yaw(j).cost < yaw(j).best.cost
            yaw(j).best.position = yaw(j).position;
            yaw(j).best.cost     = yaw(j).cost;
        end

        % Mejores globales
        if pitch(j).best.cost < gb_pitch.cost, gb_pitch = pitch(j).best; end
        if yaw(j).best.cost   < gb_yaw.cost,   gb_yaw   = yaw(j).best;   end
    end

    best_costs_pitch(it) = gb_pitch.cost;
    best_costs_yaw(it)   = gb_yaw.cost;

    if partial_log
        fprintf('\nIter %d/%d\n', it, max_iter);
        fprintf('  BestCost_pitch = %.6g | BestPos_pitch = [%g %g %g]\n', ...
            best_costs_pitch(it), gb_pitch.position(1), gb_pitch.position(2), gb_pitch.position(3));
        if ignore_yaw
            fprintf('  (Yaw ignorado en coste)  BestPos_yaw = [%g %g %g]\n', ...
                gb_yaw.position(1), gb_yaw.position(2), gb_yaw.position(3));
        else
            fprintf('  BestCost_yaw   = %.6g | BestPos_yaw   = [%g %g %g]\n', ...
                best_costs_yaw(it), gb_yaw.position(1), gb_yaw.position(2), gb_yaw.position(3));
        end
    end
end

% ===== Resultado final: publica al workspace y muestra =====
assignin('base','p1', gb_pitch.position(1));
assignin('base','p2', gb_pitch.position(2));
assignin('base','p3', gb_pitch.position(3));
assignin('base','y1', gb_yaw.position(1));
assignin('base','y2', gb_yaw.position(2));
assignin('base','y3', gb_yaw.position(3));

fprintf('\n=== RESULTADO FINAL ===\n');
fprintf('Pitch: Jp=%.6g | [p1 p2 p3]=[%g %g %g]\n', gb_pitch.cost, gb_pitch.position);
if ignore_yaw
    fprintf('Yaw:   (ignorado en coste) | [y1 y2 y3]=[%g %g %g]\n', gb_yaw.position);
else
    fprintf('Yaw:   Jy=%.6g | [y1 y2 y3]=[%g %g %g]\n', gb_yaw.cost, gb_yaw.position);
end

%% ================= FUNCIONES LOCALES =================

function v_spd = project_to_spd_vec(v)
    v = v(:);
    A = [v(1) v(2); v(2) v(3)];
    A = (A+A.')/2;
    [~,p] = chol(A);
    if p>0
        eps_shift = 1e-6;
        d = max(0, -min(eig(A))) + eps_shift;
        A = A + d*eye(2);
    end
    v_spd = [A(1,1); A(1,2); A(2,2)];
end

function [Jp, Jy] = eval_cost(mdl, vp, vy, ignore_yaw_flag, use_rmse_flag, rmse_p_name, rmse_y_name)
    % Publica p1..p3 y y1..y3
    assignin('base','p1', vp(1));
    assignin('base','p2', vp(2));
    assignin('base','p3', vp(3));

    if ~ignore_yaw_flag
        assignin('base','y1', vy(1));
        assignin('base','y2', vy(2));
        assignin('base','y3', vy(3));
    else
        % Mantén y* existentes si ignoras yaw
        try evalin('base','y1;'); catch, assignin('base','y1',1); end
        try evalin('base','y2;'); catch, assignin('base','y2',0); end
        try evalin('base','y3;'); catch, assignin('base','y3',1); end
    end

    % Fuerza releer variables del workspace
    set_param(mdl,'SimulationCommand','update');
      
    evalin('base','clear RMSE_pitch RMSE_yaw');

    % Simulación
    try
        out = sim(mdl,'ReturnWorkspaceOutputs','on');
    catch ME
        warning('Sim error: %s', ME.message);
        Jp = inf; Jy = inf; return
    end

    if use_rmse_flag
        % Lee directamente RMSE desde out o base
        rp = grab_sig(out, rmse_p_name);
        Jp = scalar_last(rp);

        if ignore_yaw_flag
            Jy = 0;
        else
            ry = grab_sig(out, rmse_y_name);
            Jy = scalar_last(ry);
        end
    else

    end
end

function e = grab_sig(out, primary_name, fallback_name)
    e = [];
    try, v = out.get(primary_name); e = unwrap(v); end
    if isempty(e) && nargin>=3
        try, v = out.get(fallback_name); e = unwrap(v); end
    end
    if isempty(e)
        try, v = evalin('base',primary_name); e = unwrap(v); end
    end
    if isempty(e) && nargin>=3
        try, v = evalin('base',fallback_name); e = unwrap(v); end
    end
end

function x = unwrap(v)
    x = [];
    if isa(v,'timeseries')
        x = double(v.Data(:));
    elseif isstruct(v)
        if isfield(v,'signals') && isfield(v.signals,'values')
            x = double(v.signals.values(:));
        elseif isfield(v,'Data')
            x = double(v.Data(:));
        end
    elseif isnumeric(v)
        x = double(v(:));
    end
end

function val = scalar_last(sig)
    % Acepta escalar o vector/serie; toma el último valor válido
    if isempty(sig) || all(isnan(sig))
        val = inf;
        return
    end
    sig = sig(:);
    idx = find(~isnan(sig) & isfinite(sig), 1, 'last');
    if isempty(idx)
        val = inf;
    else
        val = sig(idx);
    end
end

function m = mse_of(e)
    if isempty(e) || any(isnan(e))
        m = inf;
    else
        n = max(1,numel(e));
        m = sum(e.^2)/n;
    end
end
