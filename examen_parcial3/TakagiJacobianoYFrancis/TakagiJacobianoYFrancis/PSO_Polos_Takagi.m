function result = PSO_poles()
% PSO para optimizar p=[p1 p2 p3 p4] (polos deseados)
% Minimiza el RMSE de la simulación:
%   1) Publica p1..p4 en el base workspace
%   2) Ejecuta KnlKnsKz.m (debe usar p = [p1 p2 p3 p4])
%   3) Corre el modelo TakagiSugenaConControl1.slx
%   4) Lee out.RMSE (variable 'RMSE')

%% ===== CONFIG =====
mdl         = 'TakagiSugenaConControl1';   % nombre del modelo .slx sin extensión
rmse_name   = 'RMSE';                      % nombre de la señal de desempeño
gain_script = 'KnlKnsKz';                  % script que calcula K_nl,K_ns,K_z,P_*,G_*
max_iter    = 50;
pop_size    = 100;
c1 = 0.40; c2 = 0.50; w = 0.60;
partial_log = true;

% Warm-start opcional (si quieres algún p inicial conocido)
use_init   = false;           % pon true si quieres iniciar cerca de un vector dado
init_pvec = [-5 -6 -7 -8].'; % ejemplo de polos iniciales (todos en R-)
rng(1);
n = 4;                        % dimensión (p1..p4)

%% ===== SIMULINK SETUP =====
load_system(mdl);
set_param(mdl,'ReturnWorkspaceOutputs','on','SimulationMode','normal');
% set_param(mdl,'FastRestart','off');
% set_param(mdl,'SimulationCommand','update');

%% ===== SWARM STRUCT =====
empty.position      = zeros(n,1);
empty.velocity      = zeros(n,1);
empty.cost          = inf;
empty.best.position = zeros(n,1);
empty.best.cost     = inf;

swarm = repmat(empty, pop_size, 1);
gbest = empty.best; gbest.cost = inf;

%% ===== INIT SWARM =====
for i = 1:pop_size
    if use_init
        if i == 1
            p = init_pvec(:);
        else
            p = init_pvec(:) .* (1 + 0.05*randn(n,1));
        end
    else
        % Inicializa polos en el rango [-20, -1]
        p = -1 - (20-1)*rand(n,1);    % valores entre ~-1 y -20
    end

    p = project_to_stable_vec(p);
    swarm(i).position = p;
    swarm(i).velocity = zeros(n,1);
    swarm(i).cost     = eval_cost(mdl, p, rmse_name, gain_script);
    swarm(i).best.position = p;
    swarm(i).best.cost     = swarm(i).cost;

    if swarm(i).best.cost < gbest.cost
        gbest = swarm(i).best;
    end
end

%% ===== Post-init: gbest válido =====
costs = [swarm.cost];
if all(~isfinite(costs))
    error('PSO:init:AllInfCosts', ...
        'Todos los costos son inf/NaN. Verifica que "%s" produzca "%s".', mdl, rmse_name);
end
finite_idx = find(isfinite(costs), 1, 'first');
if isempty(gbest.position) || ~isfinite(gbest.cost)
    gbest = swarm(finite_idx).best;
end
gbest.position = gbest.position(:);

best_costs = nan(max_iter,1);

%% ===== LOOP PSO =====
for it = 1:max_iter
    for j = 0:pop_size-1
        idx = j+1;

        r1 = rand(n,1);
        r2 = rand(n,1);

        pj  = swarm(idx).position(:);
        vj  = swarm(idx).velocity(:);
        pbl = swarm(idx).best.position(:);
        pbg = gbest.position(:);

        % Actualiza velocidad y posición
        vj = w*vj + c1.*r1.*(pbl - pj) + c2.*r2.*(pbg - pj);
        pj = pj + vj;

        % Proyección a polos estables (reales negativos, acotados)
        pj = project_to_stable_vec(pj);
        J  = eval_cost(mdl, pj, rmse_name, gain_script);

        swarm(idx).velocity = vj;
        swarm(idx).position = pj;
        swarm(idx).cost     = J;

        if J < swarm(idx).best.cost
            swarm(idx).best.position = pj;
            swarm(idx).best.cost     = J;
        end
        if swarm(idx).best.cost < gbest.cost
            gbest = swarm(idx).best;
        end
    end

    best_costs(it) = gbest.cost;
    if partial_log
        fprintf('Iter %d/%d | BestCost=%.6g | BestPoles=[%g %g %g %g]\n', ...
            it, max_iter, best_costs(it), ...
            gbest.position(1), gbest.position(2), gbest.position(3), gbest.position(4));
    end
end

%% ===== RESULTADOS =====
p1 = gbest.position(1);
p2 = gbest.position(2);
p3 = gbest.position(3);
p4 = gbest.position(4);

assignin('base','p1',p1);
assignin('base','p2',p2);
assignin('base','p3',p3);
assignin('base','p4',p4);
assignin('base','PSO_best_cost', gbest.cost);

result.p1      = p1;
result.p2      = p2;
result.p3      = p3;
result.p4      = p4;
result.cost    = gbest.cost;
result.history = best_costs;

fprintf('\n=== RESULTADO FINAL ===\n');
fprintf('Coste=%.6g | p=[%g %g %g %g]\n', gbest.cost, p1, p2, p3, p4);
end

%% ===== FUNCIONES LOCALES =====
function v_stable = project_to_stable_vec(v)
% Asegura que los polos sean reales negativos y acotados:
%   -50 <= p_i <= -1
    v = v(:);
    % Forzar signo negativo
    v = -abs(v);
    % Cotas
    v = max(v, -50);   % no más rápidos que -50
    v = min(v, -1);    % no más cercanos a 0 que -1
    v_stable = v;
end

function J = eval_cost(mdl, pvec, rmse_name, gain_script)
    % Publica p1..p4 en base
    assignin('base','p1', pvec(1));
    assignin('base','p2', pvec(2));
    assignin('base','p3', pvec(3));
    assignin('base','p4', pvec(4));

    % Ejecuta script de ganancias (usa p=[p1 p2 p3 p4])
    try
        evalin('base', gain_script);   % corre KnlKnsKz.m en el workspace base
    catch ME
        warning('Gain script error: %s', ME.message);
        J = inf; 
        return
    end

    % Simula el modelo
    set_param(mdl,'ReturnWorkspaceOutputs','on');
    try
        out = sim(mdl,'ReturnWorkspaceOutputs','on');
    catch ME
        warning('Sim error: %s', ME.message);
        J = inf; 
        return
    end

    % Leer RMSE desde SimulationOutput o base
    e = [];
    % a) campo directo (out.RMSE)
    try
        v = out.(rmse_name); 
        e = unwrap_signal(v);
    end

    % b) método get
    if isempty(e)
        try
            v = out.get(rmse_name); 
            e = unwrap_signal(v);
        end
    end

    % c) variable en base
    if isempty(e)
        try
            v = evalin('base', rmse_name); 
            e = unwrap_signal(v);
        end
    end

    J = scalar_last(e);
end

function x = unwrap_signal(v)
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
