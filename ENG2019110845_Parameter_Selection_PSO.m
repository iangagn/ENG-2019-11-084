%% PREPROCESS

% Set verbosity
verbose = true;

% Run parameters
nb_runs = 1e3;

% Objective funtions
fun = {@sphere, @rosenbrock, @linear_step, @noisy_quartic, @foxholes};

% Domain limits
lb = [-5.12, -5.12, -5.12, -1.28, -65.536]; 
ub = abs(lb);

% Global minima
min_values = [0, 0, 0, 0, 0.998];

% Parameter settings
allparams = allcomb([0.7, 0.9, 1.2, 1.5], [0.7, 0.9, 1.3, 1.5], [0.7, 0.9, 1.3, 1.5]);
nb_paramset = size(allparams, 1);

% Initialization
nb_functions = numel(fun);
for i = 1 : nb_functions
    for k = 1  : nb_paramset
        results.(func2str(fun{i})).(['paramset', num2str(k)]).fht = nan(nb_runs, 1);
    end
end

% Parameter set generating function
fparams = @(p, lb, ub, f, minval) struct('fun',                 f,      ...
                                         'nb_dim',              2,      ...
                                         'initial_positions',   [],     ...
                                         'lower_bound',         lb,     ...
                                         'upper_bound',         ub,     ...
                                         'w',                   p(1),   ...
                                         'c1',                  p(2),   ...
                                         'c2',                  p(3),   ...
                                         'nb_particles',        20,     ...
                                         'max_iter',            2500,   ...
                                         'known_best_fitness',  minval, ...
                                         'tol',                 1e-2,   ...
                                         'positions_hist_flag', false);

%% MAIN

% Initialize progress tracker
count = 0;

% Calculate the total number of runs 
nb_total_runs = nb_functions * nb_paramset * nb_runs;

% For each objective function
for i = 1 : nb_functions
    
  % For each parameter set
  for j = 1 : nb_paramset
      
      % Calculate first hitting times
      for k = 1 : nb_runs
          count = count + 1;
          results.(func2str(fun{i})).(['paramset', num2str(j)]).fht(k) = pso(fparams(allparams(j,:), lb(i), ub(i), fun{i}, min_values(i)));
      end
      
      % Show progress
      if verbose
          disp(['Completed : ' num2str(count/nb_total_runs*100) ' %'])
      end
      
      % Calculate nanmedian (i.e. ignores NaNs)
      results.(func2str(fun{i})).(['paramset', num2str(j)]).median = nanmedian(results.(func2str(fun{i})).(['paramset', num2str(j)]).fht);
      
      % Calculate percentage of targets hit
      results.(func2str(fun{i})).(['paramset', num2str(j)]).converged_percentage = 1 - sum(isnan(results.(func2str(fun{i})).(['paramset', num2str(j)]).fht))/nb_runs;
  end
  
end

%% POSTPROCESS

for i = 1 : nb_functions-1
    
    % Create rankings structure
    rankings.(func2str(fun{i})) = nan(nb_paramset, 2);
    
    for j = 1 : nb_paramset
        for k = 1 : nb_runs
            
            % Calculate percentage of hits
            converged_percentage = (1 - sum(isnan(results.(func2str(fun{i})).(['paramset', num2str(j)]).fht))/nb_runs)*100;
            results.(func2str(fun{i})).(['paramset', num2str(j)]).converged_percentage = converged_percentage;
            
            % Calculate median and adjusted median
            median = nanmedian(results.(func2str(fun{i})).(['paramset', num2str(j)]).fht);
            results.(func2str(fun{i})).(['paramset', num2str(j)]).median = median;
            adjusted_median = median/(converged_percentage/100);
            results.(func2str(fun{i})).(['paramset', num2str(j)]).adjusted_median = adjusted_median;
           
        end
        
         % Store adjusted median for ranking
         rankings.(func2str(fun{i}))(j,1) = adjusted_median;
        
    end
    
    % Store adjusted median for ranking
    [~, rankings.(func2str(fun{i}))(:,2)] = sort(rankings.(func2str(fun{i}))(:,1),'ascend');
    
end

% Find best parameter set (i.e. minimum sum of ranks over the test set)
sum_of_ranks = zeros(nb_paramset, 1);
for i = 1 : nb_functions
    sum_of_ranks = sum_of_ranks + rankings.(func2str(fun{i}))(:,2);
end

[minrank, minrank_idx] = min(sum_of_ranks);

best_paramset = allparams(minrank_idx, :);

%% Explore results

idx = minrank_idx;
disp(['f1 = ', num2str(results.sphere.(['paramset', num2str(idx)]).adjusted_median)])
disp(['f2 = ', num2str(results.rosenbrock.(['paramset', num2str(idx)]).adjusted_median)])
disp(['f3 = ', num2str(results.linear_step.(['paramset', num2str(idx)]).adjusted_median)])
disp(['f4 = ', num2str(results.noisy_quartic.(['paramset', num2str(idx)]).adjusted_median)])
disp(['f5 = ', num2str(results.foxholes.(['paramset', num2str(idx)]).adjusted_median)])
