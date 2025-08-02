function plot_figure2(l, x, OF, sol_list)
fig = figure;
set(fig, 'Color', 'w');

ax = gca;
set(ax, 'Color', 'w');
hold on;

% plot Onsager Feynman results
plot(x, OF.^2, Color = 'Red', LineStyle = ':', LineWidth = 2)

% plot R^2 and W for each gamma
styles = {'-.', '--', ':', '-'};
for i = 1:length(sol_list)
    plot(x, sol_list(i).R.^2, 'Color', [0,0,0], 'LineStyle', styles{i}, LineWidth =  2);
end

% define strong screening limits (gamma = 0.01)
r_strong_screening = sqrt(1+4*0.01^2*l^2.*(1./x.^4));
w_strong_screening = -l^2./x.^2;

% plot R strong screening limit
plot(x, r_strong_screening.^2, Color = 'Red', LineWidth = 2')

% plot each W
for i = 1:length(sol_list)
    plot(x, sol_list(i).W, 'Color', [0.5 0.5 0.5], LineStyle = styles{i}, LineWidth = 2);
end

% plot W strong screening limit 
plot(x, w_strong_screening, Color = 'Red', LineWidth = 2, LineStyle = '--')

yline(1, '--');
yline(0);
legend({'R^2, \gamma → ∞', ...
        'R^2, \gamma = 5', ...
        'R^2, \gamma = 1', ...
        'R^2, \gamma = 0.1', ...
        'R^2, \gamma = 0.01', ...
        'Eq.(50),\gamma = 0.01',...
        'w, \gamma = 5', ...
        'w, \gamma = 1', ...
        'w, \gamma = 0.1', ...
        'w, \gamma = 0.01', ... 
        'Eq.(49)'}, ...
        'Location', 'southeast');

xlabel('\rho');
xlim([0, 7]);
ylim([-2, 1.2]);
hold off;
end
