
% Run parameters
nb_runs = 1000;

% Objective funtions
fun = {@sphere, @rosenbrock, @linear_step, @noisy_quartic, @foxholes};

% Domain limits
lb = [-5.12, -5.12, -5.12, -1.28, -65.536];
ub = abs(lb);

% Global minima
min_values = [0, 0, 0, 0, 0.998];

% Initialization
nb_functions = numel(fun);
for i = 1 : nb_functions
   results_pso_1.(func2str(fun{i})).fht = nan(nb_runs,1);
   results_pso_2.(func2str(fun{i})).fht = nan(nb_runs,1);
end

%% Run

% Parameter set generating function
fparams1 = @(lb, ub, f, minval) struct('fun',                  f,      ...
                                        'nb_dim',              2,      ...
                                        'initial_positions',   [],     ...
                                        'lower_bound',         lb,     ...
                                        'upper_bound',         ub,     ...
                                        'w',                   0.7,    ...
                                        'c1',                  0.7,    ...
                                        'c2',                  0.7,    ...
                                        'nb_particles',        20,     ...
                                        'max_iter',            5000,   ...
                                        'known_best_fitness',  minval, ...
                                        'tol',                 1e-2,   ...
                                        'positions_hist_flag', false);

fparams2 = @(lb, ub, f, minval) struct('fun',                  f,      ...
                                        'nb_dim',              2,      ...
                                        'initial_positions',   [5,-32,-32], ...
                                        'lower_bound',         lb,     ...
                                        'upper_bound',         ub,     ...
                                        'w',                   0.7,    ...
                                        'c1',                  0.7,    ...
                                        'c2',                  0.7,    ...
                                        'nb_particles',        20,     ...
                                        'max_iter',            5000,   ...
                                        'known_best_fitness',  minval, ...
                                        'tol',                 1e-2,   ...
                                        'positions_hist_flag', false);
                                                                   
count = 0;
nb_total_runs = nb_functions * nb_runs;
for i = 1 : nb_functions
    for j = 1 : nb_runs
        results_pso_1.(func2str(fun{i})).fht(j) = pso(fparams1(lb(i), ub(i), fun{i}, min_values(i)));
        results_pso_2.(func2str(fun{i})).fht(j) = pso(fparams2(lb(i), ub(i), fun{i}, min_values(i)));
        count = count + 1;
        if mod(count, 1) == 0
            disp(['Completed : ' num2str(count/nb_total_runs*100) ' %'])
        end
    end
end

%% Compare distributions

function_name = 'foxholes';

boxplot([results_pso_1.(function_name).fht, results_pso_2.(function_name).fht], {'PSO1', 'PSO2'});
set(gca, 'color', [253,245,230]/255, 'fontname', 'times', 'fontsize', 14)

text_handles = findall(gca, 'type', 'text');
for i = 1 : numel(text_handles)
   temp = get(text_handles(i)); 
   if strcmpi(temp.String, 'PSO1') || strcmpi(temp.String, 'PSO2')
       position = get(text_handles(i), 'position');
       position(2) = position(2) - 10;
       set(text_handles(i), 'fontname', 'times', 'fontsize', 12, 'position', position)
   end
end

disp(['PSO1 = ', num2str(sum(isnan(results_pso_1.(function_name).fht)))]); 
disp(['PSO2 = ', num2str(sum(isnan(results_pso_2.(function_name).fht)))]); 

nanmedian(results_pso_1.(function_name).fht)
nanmedian(results_pso_2.(function_name).fht)

%% Test

pso1 = results_pso_1.(function_name).fht(~isnan(results_pso_1.(function_name).fht));
pso2 = results_pso_2.(function_name).fht(~isnan(results_pso_2.(function_name).fht));

[h_ks, p_ks] = kstest2(pso1, pso2, 'Alpha', 0.05);
[p_wmw, h_wmw] = ranksum(pso1, pso2, 'Alpha', 0.05);

