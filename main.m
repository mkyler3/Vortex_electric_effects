 clear; clc;

% set parameters
gamma_list = [5, 1, 0.1, 0.01]; % screening strengths
l = 1; 
xmax = 100; % system length
steps = 100000; % resolution

% solve ODEs for each gamma
sol_list = struct();
for i = 1:length(gamma_list)
    gamma = gamma_list(i);
    [x, R_square, W] = Utils.electric_effects(l,gamma, steps, xmax);
    sol_list(i).gamma = gamma;
    sol_list(i).x = x;
    sol_list(i).R = R_square;
    sol_list(i).W = W;
end

% get Onsager-Feynman (gamma -> inf) solution
[x, OF] = Utils.onsager_feynman(l, steps, xmax);  

% plot Figure 2
Utils.plot_figure2(l,x, OF, sol_list);