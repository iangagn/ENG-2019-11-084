function output = bat_roulette_wheel(input)
% BAT.m is a parallel implementation of the Bat algorithm as described in
% [Yang, 2010] with the exception that it uses roulette wheel selection.
%
% NOTE: The gamma parameter is no longer used.
%
% Author: Iannick Gagnon (iannick.gagnon.1@ens.etsmtl.ca)
%
%% Input checks
expected_fieldnames = {     ...
    'fun';                  ... % The objective function
    'nb_dim';               ... % The number of dimensions of the hypercube domain
    'initial_positions';    ... % The initial positions [OPTIONAL and random uniform by default]
    'lower_bound';          ... % Lower bound of the hypercube domain
    'upper_bound';          ... % Upper boind of the hypercube domain
    'fmin';                 ... % Min frequency value in BA (equation 2 in (Yang, 2010))
    'fmax';                 ... % Max frequency value in BA (equation 2 in (Yang, 2010))
    'r0';                   ... % Initial pulse rate in BA (equation 6 in (Yang, 2010))
    'a';                    ... % Alpha parameter in BA (equation 6 in (Yang, 2010))
    'g';                    ... % Gamma parameter in BA (equation 6 in (Yang, 2010))
    'loudness';             ... % Loudness parameter in BA (equation 5 in (Yang, 2010))
    'epsilon';              ... % Epsilon in BA (equation 5 in (Yang, 2010))
    'nb_bats';              ... % Number of agents/bats
    'max_iter';             ... % Maximum number of iterations
    'known_best_fitness';   ... % Known best fitness value of objective function
    'tol';                  ... % Tolerance such that (Value - tol == Value) = True
    'positions_hist_flag'};     % Flag that triggers the tracking of the positions           

% Make sure that the INPUT variable contains all of the required parameters
input_fieldnames = fieldnames(input);
assert(all(strcmpi(expected_fieldnames, input_fieldnames)), 'ERROR: Missing input(s)');

% Transform structure fields-values into workspace variables
for i = 1 : numel(expected_fieldnames)
    evalc([input_fieldnames{i} ' = input.' input_fieldnames{i}]);
end

% Make sure the objective function is a function_handle object
assert(isa(fun, 'function_handle'), 'ERROR: Objective function must be of function_handle type')

% Make sure the number of dimensions is greater than zero
assert(nb_dim > 0, 'ERROR: Objective function dimensions must be > 0');

% Make sure the lower bound is less than the upper bound
assert(lower_bound < upper_bound, 'ERROR: The lower bound cannot be greater than the upper bound')

%% Initialization

% Assign initial positions (1) randomly or (2) equal to input values
if isempty(initial_positions)
    
    % Assign random uniform initial positions
    positions = lower_bound * ones(nb_bats, nb_dim) + rand(nb_bats, nb_dim) * (upper_bound-lower_bound);
 
else
     
    l = initial_positions(1);
    gg = [initial_positions(2), initial_positions(3)];
    
    positions = lower_bound * ones(nb_bats, nb_dim) + rand(nb_bats, nb_dim) * (upper_bound-lower_bound);
    
    for i = 1 : nb_dim
        
        min1 = gg(i) - l;
        max1 = gg(i) + l;
        
        min2 = gg(2) - l;
        max2 = gg(2) + l;
        
        bad_idx = find((~(double(positions(:,1) < min1) + double(positions(:,1) > max1)) + ~(double(positions(:,2) < min2) + double(positions(:,2) > max2)))==2);
        nb_bad = numel(bad_idx);
        
        while ~isempty(bad_idx)
            
            positions(bad_idx, :) = lower_bound * ones(nb_bad, 2) + rand(nb_bad, 2) * (upper_bound-lower_bound);
            bad_idx = find((~(double(positions(:,1) < min1) + double(positions(:,1) > max1)) + ~(double(positions(:,2) < min2) + double(positions(:,2) > max2)))==2);
            nb_bad = numel(bad_idx);
            
        end
        
    end
    
    %     % Make sure the initial positions array dimensions are coherent with the number of agents and objective function dimensions
    %     assert(numel(initial_positions) == 2*nb_dim, 'ERROR: Invalid input initial positions bounds');
    %
    %     positions =  ones(nb_bats, nb_dim);
    %     for i = 1 : nb_dim
    %         positions(:,i) = initial_positions(2*i-1) * positions(:,i) + rand(nb_bats, 1) * (initial_positions(2*i)-initial_positions(2*i-1));
    %     end
    
end

% Initialize fitness values
fitnesses = fun(positions);

% Initialize number of objective function evaluations
nb_fun_eval = nb_bats;

% Initialize global best fitness value
global_best_fitness = min(fitnesses);

% Initialize global best positions
global_best_position = positions(find(fitnesses == min(fitnesses), 1, 'first'),:);

% Initialize first hitting time
first_hitting_time = NaN;

% Initialize loudnesses
A = loudness * ones(nb_bats, 1);

% Initialize velocities
velocities = zeros(nb_bats, nb_dim);

% Initialize number of epochs
nb_epochs = 0;

% Initialize positions history
positions_history = {};

%% Main loop

for iteration = 1 : max_iter
    
    % Increment number of epochs
    nb_epochs = nb_epochs + 1;
    
    % Random walk around global best using roulette wheel selection
    random_walk_indexes = find(chance([0.5, 0.5], nb_bats)-1);
    if ~isempty(random_walk_indexes)
        nb_random_walks = numel(random_walk_indexes);
        positions(random_walk_indexes, :) = repmat(global_best_position, nb_random_walks, 1) + epsilon * randn(nb_random_walks, nb_dim);
    end
    
    
    % New positions
    frequency = fmin + (fmax - fmin) * rand(nb_bats, nb_dim);
    dist_with_best = repmat(global_best_position, nb_bats, 1) - positions ;
    velocities = velocities + dist_with_best .* frequency;
    positions_new = positions + velocities;
    
    % Random uniform reinitialization of agents that are out of bounds
    [rlb, clb] = find(positions_new < lower_bound);
    [rub, cub] = find(positions_new > upper_bound);
    row = [rlb; rub];
    col = [clb; cub];
    if ~isempty(row)
        nb_out_of_bounds = numel(row);
        for i = 1 : nb_out_of_bounds
            positions_new(row(i),col(i)) = lower_bound + rand() * (upper_bound-lower_bound);
        end
    end
    
    % Calculate new fitnesses
    fitnesses_new = fun(positions_new);
    [new_best, best_index] = min(fitnesses_new);
    
    % Increment number of objetive function evaluations
    nb_fun_eval = nb_fun_eval + nb_bats;
    
    % Update if improving AND not-too-loud
    if min(fitnesses_new) < global_best_fitness && rand() < A(best_index)
        
        % Accept new positions
        positions = positions_new;
        
        % Update loudness and pulse rate
        not_too_loud_indexes = find(rand(nb_bats, 1) < A);
        A(not_too_loud_indexes) = a * A(not_too_loud_indexes);
 
        % Update global best
        global_best_fitness = new_best;
        global_best_position = positions(best_index,:);
        
    end
    
    % Update positions history
    if positions_hist_flag
        positions_history{end + 1} = positions; %#ok
    end
    
    % Check for convergence
    if abs(global_best_fitness - known_best_fitness) < tol
        first_hitting_time = nb_fun_eval;
        break
    end

end

%% Output

output = first_hitting_time;

end