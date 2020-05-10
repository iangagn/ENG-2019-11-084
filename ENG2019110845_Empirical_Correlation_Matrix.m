%% Calculate correlation coefficients

functions = {@sphere, @rosenbrock, @linear_step, @noisy_quartic, @foxholes};
nb_functions = numel(functions);

% 
corr_mat = empirical_correlation(functions, {[-5.12, 5.12], [-5.12, 5.12], [-5.12, 5.12], [-1.28, 1.28], [-65.536, 65.536]}, 1000);

%% Plot

% Image + colorbar
h = imagesc(corr_mat);
set(h, 'alphadata', 0.6)
colormap(jet);

% Labels
labels = {'$f_{1}$','$f_{2}$','$f_{3}$', '$f_{4}$', '$f_{5}$'};
set(gca,'xtick', 1:nb_functions,'xticklabel',labels, 'ytick', 1:nb_functions, 'yticklabel',labels, 'fontname', 'times', 'fontsize', 18, 'TickLabelInterpreter', 'latex'); 

% Show values as text
[rows,cols] = size(corr_mat);
for i = 1 : rows
    for j = 1 : cols
        textHandles(j,i) = text(j,i,num2str(corr_mat(i,j)),'horizontalAlignment','center', 'fontname', 'times', 'fontsize', 18); %#ok
    end
end

% Calculate mean correlation coefficient
lower_triangular = tril(corr_mat, -1);
mean_corr = sum(sum(lower_triangular)) / sum(1:nb_functions-1);

