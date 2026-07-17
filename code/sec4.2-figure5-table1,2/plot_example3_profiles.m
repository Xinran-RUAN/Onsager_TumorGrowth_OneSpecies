% plot_example3_profiles.m
% Plot final density and pressure profiles for Example 3.

clear; clc; close all;

% Parameters
outdir = 'example3';
gammaList = [3, 20, 200];
lineWidth = 2.7;
lineColors = {
   
        [0.47, 0.67, 0.19]   % green
    [0.93, 0.69, 0.13]   % yellow
    [0.85, 0.10, 0.10]   % red
 [0.00, 0.45, 0.74]   % blue
};
lineStyles = {'--', '-.', ':'};   % finite gamma profiles
hsStyle = '-';                    % Hele-Shaw profile
outputFile = 'Figure8_example3_profiles.png';

data = cell(length(gammaList), 1);

for ig = 1:length(gammaList)
    gamma = gammaList(ig);
    files = dir(fullfile(outdir, sprintf('example3_gamma%d_N*.mat', gamma)));
    if isempty(files)
        error('No data file found for gamma = %d.', gamma);
    end
    [~, idx] = max([files.datenum]);
    data{ig} = load(fullfile(outdir, files(idx).name));
end

% Use the last data file to define the Hele-Shaw limit at final time
S = data{end};
x  = S.x;
T  = S.T;
R0 = S.R0;

R_T = asinh(exp(T) * sinh(R0));
p_hs = 1 - cosh(x) / cosh(R_T);
p_hs = max(p_hs, 0);
n_hs = double(abs(x) <= R_T);

fig = figure;
set(fig, 'Color', 'w');
set(fig, 'Position', [100, 100, 1100, 420]);

legendText = cell(length(gammaList) + 1, 1);
for ig = 1:length(gammaList)
    legendText{ig} = ['$\gamma=', num2str(gammaList(ig)), '$'];
end
legendText{end} = 'Hele--Shaw';

% =====================================================
% Density
% =====================================================
ax1 = subplot(1,2,1);
hold on;

h1 = gobjects(length(gammaList) + 1, 1);

h1(end) = plot(x, n_hs, ...
    'Color', lineColors{length(gammaList) + 1}, ...
    'LineStyle', hsStyle, ...
    'LineWidth', 2.7);

for ig = 1:length(gammaList)
    S = data{ig};
    h1(ig) = plot(S.x, S.sol_n(:,end), ...
        'Color', lineColors{ig}, ...
        'LineStyle', lineStyles{ig}, ...
        'LineWidth', lineWidth);
end

xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$n$', 'Interpreter', 'latex', 'FontSize', 22);
title(['Density at $T=', num2str(T, '%.0f'), '$'], ...
    'Interpreter', 'latex', 'FontSize', 22);

set(gca, ...
    'FontSize', 22, ...
    'LineWidth', 1.5, ...
    'TickLabelInterpreter', 'latex');

box on;
grid on;

xlim([-5, 5]);
xticks([-5, -2.5, 0, 2.5, 5]);
set(gca, 'XTickLabel', {'$-5$', '$-2.5$', '$0$', '$2.5$', '$5$'});

ylim([0, 1.2]);
yticks([0:0.3:1.2]);
% set(gca, 'YTickLabel', {'$0$', '$0.2$', '$0.4$', '$0.6$', '$0.8$', '1'});

lgd1 = legend(h1, legendText, ...
    'Interpreter', 'latex', ...
    'Location', 'northeast', ...
    'Box', 'on', ...
    'FontSize', 18);
set(lgd1, 'LineWidth', 1.5, 'EdgeColor', 'k', 'Color', 'w');

% =====================================================
% Pressure
% =====================================================
ax2 = subplot(1,2,2);
hold on;

h2 = gobjects(length(gammaList) + 1, 1);

h2(end) = plot(x, p_hs, ...
    'Color', lineColors{length(gammaList) + 1}, ...
    'LineStyle', hsStyle, ...
    'LineWidth', 2.7);

for ig = 1:length(gammaList)
    S = data{ig};
    h2(ig) = plot(S.x, S.sol_p(:,end), ...
        'Color', lineColors{ig}, ...
        'LineStyle', lineStyles{ig}, ...
        'LineWidth', lineWidth);
end

xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$p$', 'Interpreter', 'latex', 'FontSize', 22);
title(['Pressure at $T=', num2str(T, '%.0f'), '$'], ...
    'Interpreter', 'latex', 'FontSize', 22);

set(gca, ...
    'FontSize', 22, ...
    'LineWidth', 1.5, ...
    'TickLabelInterpreter', 'latex');

box on;
grid on;

xlim([-5, 5]);
xticks([-5, -2.5, 0, 2.5, 5]);
set(gca, 'XTickLabel', {'$-5$', '$-2.5$', '$0$', '$2.5$', '$5$'});

ylim([0, 1.2]);
yticks([0:0.3:1.2]);
% set(gca, 'YTickLabel', {'$0$', '$0.2$', '$0.4$', '$0.6$', '$0.8$', '1'});

lgd2 = legend(h2, legendText, ...
    'Interpreter', 'latex', ...
    'Location', 'northeast', ...
    'Box', 'on', ...
    'FontSize', 18);
set(lgd2, 'LineWidth', 1.5, 'EdgeColor', 'k', 'Color', 'w');

% subplot positions
set(ax1, 'Position', [0.06, 0.16, 0.41, 0.72]);
set(ax2, 'Position', [0.57, 0.16, 0.41, 0.72]);

print(fig, fullfile(outdir, outputFile), '-dpng', '-r300');

