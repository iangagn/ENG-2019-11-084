
% Run parameters
nb_runs = 10;

% Objective funtions
fun = {@happycat};

% Domain limits
lb = -6.4;
ub = 6.35;

% Global minimum
min_values = 0;

% Initialization
nb_functions = numel(fun);
for i = 1 : nb_functions
   results1.(func2str(fun{i})).fht = nan(nb_runs,1);
   results2.(func2str(fun{i})).fht = nan(nb_runs,1);
end

%% Run

% Parameter set generating function
fparams1 = @(lb, ub, f, minval) struct('fun',                  f,      ...
                                       'nb_dim',               2,      ...
                                       'initial_positions',    [],     ...
                                       'lower_bound',          lb,     ...
                                       'upper_bound',          ub,     ...
                                       'fmin',                 0,      ...
                                       'fmax',                 2,      ...
                                       'r0',                   0.7,    ...
                                       'a',                    0.9,    ...
                                       'g',                    0.1,    ...
                                       'loudness',             1,      ...
                                       'epsilon',              1e-3,   ...
                                       'nb_bats',              50,     ...
                                       'max_iter',             10000,   ...
                                       'known_best_fitness',   minval, ...
                                       'tol',                  1e-1,   ...
                                       'positions_hist_flag',  true);           

fparams2 = @(lb, ub, f, minval) struct( 'fun',                  f,      ...
                                        'nb_dim',              2,      ...
                                        'initial_positions',   [],     ...
                                        'lower_bound',         lb,     ...
                                        'upper_bound',         ub,     ...
                                        'w',                   0.7,    ...
                                        'c1',                  0.7,    ...
                                        'c2',                  0.7,    ...
                                        'nb_particles',        50,     ...
                                        'max_iter',            10000,   ...
                                        'known_best_fitness',  minval, ...
                                        'tol',                 1e-1,   ...
                                        'positions_hist_flag', true); 
                                   
count = 0;
nb_total_runs = nb_functions * nb_runs;
for i = 1 : nb_functions
    for j = 1 : nb_runs
        results1.(func2str(fun{i})).fht(j) = bat(fparams1(lb(i), ub(i), fun{i}, min_values(i)));
        results2.(func2str(fun{i})).fht(j) = pso(fparams2(lb(i), ub(i), fun{i}, min_values(i)));
        count = count + 1;
        if mod(count, 1) == 0
            disp(['Completed : ' num2str(count/nb_total_runs*100) ' %'])
        end
    end
end

%% Compare distributions

function_name = 'happycat';

boxplot([results1.(function_name).fht, results2.(function_name).fht], {'BA', 'PSO'});
set(gca, 'color', [253,245,230]/255, 'fontname', 'times', 'fontsize', 14)

text_handles = findall(gca, 'type', 'text');
for i = 1 : numel(text_handles)
   temp = get(text_handles(i)); 
   if strcmpi(temp.String, 'BA') || strcmpi(temp.String, 'PSO')
       position = get(text_handles(i), 'position');
       position(2) = position(2) - 10;
       set(text_handles(i), 'fontname', 'times', 'fontsize', 12, 'position', position)
   end
end

disp(['BA = ', num2str(sum(isnan(results1.(function_name).fht)))]); 
disp(['PSO = ', num2str(sum(isnan(results2.(function_name).fht)))]); 

nanmedian(results1.(function_name).fht)
nanmedian(results2.(function_name).fht)

%% Test

ba = results1.(function_name).fht(~isnan(results1.(function_name).fht));
pso = results2.(function_name).fht(~isnan(results2.(function_name).fht));

[h, p] = kstest2(ba, pso);
