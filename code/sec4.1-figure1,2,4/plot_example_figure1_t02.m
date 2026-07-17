function plot_example_figure1_t02()
% plot_example_figure1(plotT)
%
% Plot Barenblatt profiles for gamma = 3 and gamma = 20.
% If plotT is omitted, the default plotting time below is used.
% Use plotT = [] to plot the final saved time.

clc; close all;

outdir = 'example1_t02';

% Plot time. Use [] to plot the final saved time.
% For example, set plotT = 0.05 to plot the profiles at t = 0.05.
plotT = 1;

S1 = load(fullfile(outdir, 'example1_gamma3.mat'));
S2 = load(fullfile(outdir, 'example1_gamma20.mat'));

fig = figure;
set(fig, 'Color', 'w');
set(fig, 'Position', [100, 100, 1100, 420]);

% -------- left: gamma = 3 --------
ax1 = subplot(1, 2, 1);
plot_one_case(S1, plotT);

% -------- right: gamma = 20 --------
ax2 = subplot(1, 2, 2);
plot_one_case(S2, plotT);

% manually adjust subplot positions
set(ax1, 'Position', [0.06, 0.16, 0.40, 0.74]);
set(ax2, 'Position', [0.57, 0.16, 0.40, 0.74]);

print(fig, fullfile(outdir, 'Figure1_example1_profiles.png'), '-dpng', '-r300');


function plot_one_case(S, plotT)

x = S.x;
gamma = S.gamma;
[timeIndex, T] = pick_plot_time(S, plotT);

n_num = S.sol_n(:, timeIndex);
n_ex  = barenblatt_exact(x, T, gamma, S.t0, S.Cbar);

% ----- draw exact solution first -----
h_ex = plot(x, n_ex, 'r-', ...
    'LineWidth', 2.7);
hold on;

% ----- draw numerical solution second -----
h_num = plot(x, n_num, 'bo', ...
    'MarkerSize', 7, ...
    'LineWidth', 1.5, ...
    'MarkerFaceColor', 'none');

xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$n$', 'Interpreter', 'latex', 'FontSize', 22);

title(['$\gamma = ', num2str(gamma), ',\quad T = ', num2str(T), '$'], ...
    'Interpreter', 'latex', 'FontSize', 22);

lgd = legend([h_num, h_ex], ...
    {'$\mathrm{Numerical}$', '$\mathrm{Barenblatt}$'}, ...
    'Interpreter', 'latex', ...
    'Location', 'northeast', ...
    'Box', 'on', ...
    'FontSize', 22);

set(lgd, ...
    'LineWidth', 1.5, ...
    'EdgeColor', 'k', ...
    'Color', 'w');

box on;
grid on;

set(gca, ...
    'FontSize', 22, ...
    'LineWidth', 1.5, ...
    'TickLabelInterpreter', 'latex');

xlim([-10, 10]);
xticks([-5, -2.5, 0, 2.5, 5]);

if gamma == 3
    ylim([0, 1]);
    yticks([0:0.2:1.0]);
elseif gamma == 20
    ylim([0, 1]);
    yticks([0:0.2:1.0]);
else
    ymax = max([n_num(:); n_ex(:)]);
    ymax_plot = 1.10 * ymax;
    ylim([0, ymax_plot]);
    yticks(nice_ticks(0, ymax_plot, 5));
end

end


function [timeIndex, Tplot] = pick_plot_time(S, plotT)

nt = size(S.sol_n, 2) - 1;

if isfield(S, 'time') && ~isempty(S.time)
    timeVec = S.time(:);
elseif isfield(S, 'tau') && ~isempty(S.tau)
    timeVec = (0:nt)' * S.tau;
elseif isfield(S, 'T') && ~isempty(S.T)
    timeVec = linspace(0, S.T, nt + 1).';
else
    error('Cannot determine saved time levels from the data file.');
end

if numel(timeVec) ~= nt + 1
    error('The time vector length is inconsistent with sol_n.');
end

if isempty(plotT)
    targetT = timeVec(end);
else
    targetT = plotT;
end

tol = 100 * eps(max(1, max(abs(timeVec))));
if targetT < timeVec(1) - tol || targetT > timeVec(end) + tol
    error('Requested plotT = %.16g is outside the saved time interval [%.16g, %.16g].', ...
        targetT, timeVec(1), timeVec(end));
end

[dtMin, timeIndex] = min(abs(timeVec - targetT));
Tplot = timeVec(timeIndex);

if ~isempty(plotT) && dtMin > max(tol, 1e-12)
    fprintf('Requested plotT = %.8g; using closest saved time T = %.8g.\n', ...
        targetT, Tplot);
end

end


function ticks = nice_ticks(a, b, ntarget)

if b <= a
    ticks = a;
    return;
end

raw_step = (b - a) / max(ntarget - 1, 1);
pow10 = 10^floor(log10(raw_step));
r = raw_step / pow10;

if r < 1.5
    step = 1 * pow10;
elseif r < 3
    step = 2 * pow10;
elseif r < 7
    step = 5 * pow10;
else
    step = 10 * pow10;
end

tick_min = ceil(a / step) * step;
tick_max = floor(b / step) * step;

ticks = tick_min:step:tick_max;

if numel(ticks) < 3
    ticks = linspace(a, b, ntarget);
end

end
end