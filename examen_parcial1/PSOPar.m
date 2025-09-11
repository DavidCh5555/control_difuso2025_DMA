function [Kp_best, Kd_best, gbest_cost, history] = pso_fuzzy_kp_kd()

    num_particulas   = 5;
    num_componentes  = 2;   
    max_iter         = 5;
    c1               = 1.6;
    c2               = 1.9;
    w                = 1; 


    Kp_min = 0.0;  Kp_max = 50.0;
    Kd_min = 0.0;  Kd_max = 10.0;


    vmax = [ (Kp_max-Kp_min)*0.25, (Kd_max-Kd_min)*0.25 ];

    init_values = zeros(num_particulas, num_componentes);

    rng('shuffle');
    X = init_values;                
    V = zeros(num_particulas,2);    
    if all(V(:)==0)
        V = 0.01*(rand(num_particulas,2)-0.5);
    end

    X(:,1) = min(max(X(:,1), Kp_min), Kp_max);
    X(:,2) = min(max(X(:,2), Kd_min), Kd_max);

    J = inf(num_particulas,1);
    for i = 1:num_particulas
        J(i) = cost_function_sim(X(i,:));
    end

    pbest_pos  = X;
    pbest_cost = J;

    [gbest_cost, idx] = min(J);
    gbest_pos = X(idx,:);

    history = zeros(max_iter,1);

    for it = 1:max_iter
        for i = 1:num_particulas
            r1 = rand; r2 = rand;

            V(i,:) = w*V(i,:) ...
                   + c1*r1*(pbest_pos(i,:) - X(i,:)) ...
                   + c2*r2*(gbest_pos      - X(i,:));

            X(i,:) = X(i,:) + V(i,:);

            Ji = cost_function_sim(X(i,:));

            if Ji < pbest_cost(i)
                pbest_cost(i) = Ji;
                pbest_pos(i,:) = X(i,:);
            end

            if Ji < gbest_cost
                gbest_cost = Ji;
                gbest_pos  = X(i,:);
            end
        end
        history(it) = gbest_cost;

    end
    Kp_best = gbest_pos(1);
    Kd_best = gbest_pos(2);

    assignin('base','Kp_best',Kp_best);
    assignin('base','Kd_best',Kd_best);
    assignin('base','gbest_cost',gbest_cost);
    assignin('base','history',history);

end

function e_x = cost_function_sim(x)
    assignin('base','Kp',x(1));
    assignin('base','Kd',x(2));
    errorHandler = [];
    try
        out = sim('PSOFuzzyPar.slx');
        track_error = out.track_error;
    catch e
        if isa(e,'MSLException')
            disp('Error detected on simulation');
            try
                errorHandler = e.handles{1};
            catch
                errorHandler = [];
            end
            if isempty(errorHandler)
                disp('Terminating simulation');
                error('Terminated correctly');
            end
            track_error = inf;
        else
            track_error = inf;
        end
    end

    if isa(track_error,'timeseries')
        track_error = track_error.Data;
    elseif isstruct(track_error)
        if isfield(track_error,'signals') && isfield(track_error.signals,'values')
            track_error = track_error.signals.values;
        elseif isfield(track_error,'Data')
            track_error = track_error.Data;
        end
    end
    if isempty(track_error) || any(isnan(track_error(:)))
        e_x = inf;
        return;
    end
    if isscalar(track_error) && isinf(track_error)
        e_x = inf;
        return;
    end
    if size(track_error,2) > 1
        track_error = track_error(:,1);
    end
    N = size(track_error,1);
    e_x = sum(track_error.^2) / max(N,1);
end