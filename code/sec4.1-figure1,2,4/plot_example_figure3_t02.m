clear; clc; close all;

outdir = 'example1_t02';

blueLineWidth = 2.5;
redLineWidth = 2.5;
energyPlotMaxPoints = 800;
savePdfFigure = true;

S1 = load(fullfile(outdir, 'example1_gamma3.mat'));
S2 = load(fullfile(outdir, 'example1_gamma20.mat'));

mass_rel_1 = (S1.mass - S1.mass(1)) / S1.mass(1);
mass_rel_2 = (S2.mass - S2.mass(1)) / S2.mass(1);

energy_rel_1 = S1.energy / S1.energy(1);
energy_rel_2 = S2.energy / S2.energy(1);
idxEnergy1 = select_plot_indices(S1.time, energyPlotMaxPoints);
idxEnergy2 = select_plot_indices(S2.time, energyPlotMaxPoints);

fig = figure;
set(fig, 'Color', 'w');
set(fig, 'Position', [100, 100, 1100, 420]);
set(fig, 'Renderer', 'painters');

% =====================================================
% left: relative mass change
% =====================================================
ax1 = subplot(1,2,1);

hMass1 = plot(S1.time, mass_rel_1, 'b-', 'LineWidth', blueLineWidth); hold on;
hMass2 = plot(S2.time, mass_rel_2, 'r--', 'LineWidth', redLineWidth);

xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('relative mass variation', 'Interpreter', 'latex', 'FontSize', 20);
title('Mass conservation', 'Interpreter', 'latex', 'FontSize', 20);

lgd1 = legend([hMass1, hMass2], ...
    {'$\gamma = 3$', '$\gamma = 20$'}, ...
    'Interpreter', 'latex', ...
    'Location', 'northeast', ...
    'Box', 'on', ...
    'FontSize', 20);
set(lgd1, 'LineWidth', 1.8, 'EdgeColor', 'k', 'Color', 'w');

set(gca, ...
    'FontSize', 22, ...
    'LineWidth', 1.5, ...
    'TickLabelInterpreter', 'latex');

box on;
grid on;

xlim([0, 100]);
set(gca, 'XTick', [0: 20: 100]);
set(gca, 'XTickLabel', {'$0$', '$20$', '$40$', '$60$', '$80$', '$100$'});
% 
 ylim([-1.2e-15, 1.e-15]);
set(gca, 'YTick', -1e-15:5e-16:1e-15);
%set(gca, 'YTickLabel', ...
 %   {'$0$', '$0.2$', '$0.4$', '$0.6$', '$0.8$', '$1$'});


% Axis range is left close to the original MATLAB figure.


% =====================================================
% right: energy decay
% =====================================================
ax2 = subplot(1,2,2);

hEnergy1 = plot(S1.time(idxEnergy1), energy_rel_1(idxEnergy1), ...
    'b-', 'LineWidth', blueLineWidth); hold on;
hEnergy2 = plot(S2.time(idxEnergy2), energy_rel_2(idxEnergy2), ...
    'r--', 'LineWidth', redLineWidth);

xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$E_h(t)/E_h(0)$', 'Interpreter', 'latex', 'FontSize', 20);
title('Energy dissipation', 'Interpreter', 'latex', 'FontSize', 20);

lgd2 = legend([hEnergy1, hEnergy2], ...
    {'$\gamma = 3$', '$\gamma = 20$'}, ...
    'Interpreter', 'latex', ...
    'Location', 'northeast', ...
    'Box', 'on', ...
    'FontSize', 20);
set(lgd2, 'LineWidth', 1.8, 'EdgeColor', 'k', 'Color', 'w');

set(gca, ...
    'FontSize', 22, ...
    'LineWidth', 1.5, ...
    'TickLabelInterpreter', 'latex');

box on;
grid on;

xlim([0, 100]);
set(gca, 'XTick', [0:20:100]);
set(gca, 'XTickLabel', {'$0$', '$20$', '$40$', '$60$', '$80$', '$100$'});

ylim([0., 1.0]);
set(gca, 'YTick', 0.:0.2:1.0);
set(gca, 'YTickLabel', ...
    {'$0$', '$0.2$', '$0.4$', '$0.6$', '$0.8$', '$1$'});

% subplot positions
set(ax1, 'Position', [0.07, 0.16, 0.38, 0.74]);
set(ax2, 'Position', [0.59, 0.16, 0.38, 0.74]);

print(fig, fullfile(outdir, 'Figure3_example1_mass_energy.png'), '-dpng', '-r300');

if savePdfFigure
    exportgraphics(fig, fullfile(outdir, 'Figure3_example1_mass_energy.pdf'), ...
        'ContentType', 'vector');
end

function idx = select_plot_indices(x, maxPoints)

    n = numel(x);

    if isempty(maxPoints) || maxPoints <= 0 || n <= maxPoints
        idx = 1:n;
        return;
    end

    idx = unique(round(linspace(1, n, maxPoints)));
end
