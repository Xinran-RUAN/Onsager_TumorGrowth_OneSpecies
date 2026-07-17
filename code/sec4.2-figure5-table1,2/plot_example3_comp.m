% plot_example3_complementarity.m
% Plot the final-time complementarity defect for Example 3.
%
% Comp_gamma(T) = h * sum_i | p_i^K * (1 - n_i^K) |.
%
% This script uses the same data files as plot_example3_errors.m.

clear; clc; close all;

outdir = 'example3';

gammaList = [10, 20, 40, 80, 160, 320];
Nfixed = 6400;

showReferenceLine = true;
saveFigure = true;

comp_final = zeros(length(gammaList), 1);
TuseList = zeros(length(gammaList), 1);

for ig = 1:length(gammaList)

    gamma = gammaList(ig);

    filePattern = sprintf('example3_gamma%d_N%d*.mat', gamma, Nfixed);
    files = dir(fullfile(outdir, filePattern));

    if isempty(files)
        error('No data file found for gamma = %d and N = %d.', gamma, Nfixed);
    end

    % If there are multiple files, use the newest one.
    [~, idx] = max([files.datenum]);
    matfileName = fullfile(outdir, files(idx).name);

    fprintf('Reading gamma = %d, N = %d from %s\n', ...
        gamma, Nfixed, files(idx).name);

    S = load(matfileName);

    if ~isfield(S, 'sol_n') || ~isfield(S, 'sol_p')
        error('File %s must contain sol_n and sol_p.', files(idx).name);
    end

    % Read mesh size.
    if isfield(S, 'h')
        h = S.h;
    elseif isfield(S, 'x')
        x = S.x;
        x = x(:);
        h = mean(diff(x));
    else
        h = 10 / Nfixed;  % fallback for the current [-5,5] test
    end

    % Read final-time numerical density and pressure.
    sz_n = size(S.sol_n);
    sz_p = size(S.sol_p);

    if sz_n(1) == Nfixed
        n_final = S.sol_n(:, sz_n(2));
    elseif sz_n(2) == Nfixed
        n_final = S.sol_n(sz_n(1), :).';
    else
        error('Unexpected size of sol_n in %s.', files(idx).name);
    end

    if sz_p(1) == Nfixed
        p_final = S.sol_p(:, sz_p(2));
    elseif sz_p(2) == Nfixed
        p_final = S.sol_p(sz_p(1), :).';
    else
        error('Unexpected size of sol_p in %s.', files(idx).name);
    end

    n_final = n_final(:);
    p_final = p_final(:);

    % Final time, only for printing and title.
    if isfield(S, 'time')
        time = S.time;
        Tuse = time(end);
    elseif isfield(S, 'T')
        Tuse = S.T;
    else
        Tuse = NaN;
    end
    TuseList(ig) = Tuse;

    % Complementarity defect.
    comp_final(ig) = h * sum(abs(p_final .* (1 - n_final)));

    fprintf('  T = %.8e, Comp = %.8e\n', Tuse, comp_final(ig));
end

fprintf('\nSummary\n');
fprintf(' gamma        T              Comp_gamma(T)\n');
for ig = 1:length(gammaList)
    fprintf('%6d   %12.6e   %12.6e\n', ...
        gammaList(ig), TuseList(ig), comp_final(ig));
end

% ========================================================================
% Build a reference line with slope -1, only as a visual guide.
% ========================================================================

if showReferenceLine
    refSlope = -1;
    iRef = 2;  % anchor at gamma = 20
    gamma_ref = gammaList(iRef);
    y_ref = 0.7 * comp_final(iRef);
    refLine = y_ref * (gammaList / gamma_ref).^refSlope;
end

% ========================================================================
% Plot
% ========================================================================

fig = figure;
set(fig, 'Color', 'w');
set(fig, 'Position', [100, 100, 700, 480]);

loglog(gammaList, comp_final, 'o-', ...
    'LineWidth', 2.2, ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', [0.00, 0.45, 0.74], ...
    'Color', [0.00, 0.45, 0.74]);
hold on;

if showReferenceLine
    loglog(gammaList, refLine, 'k--', ...
        'LineWidth', 2.0);
end

grid on;
grid minor;
box on;

xlabel('$\gamma$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$\mathrm{Comp}_{\gamma}(T)$', ...
    'Interpreter', 'latex', 'FontSize', 18);

if ~isnan(TuseList(end))
    title(['$T=', num2str(TuseList(end), '%.3g'), '$'], ...
        'Interpreter', 'latex', 'FontSize', 18);
end

set(gca, ...
    'FontSize', 18, ...
    'LineWidth', 1.0, ...
    'XScale', 'log', ...
    'YScale', 'log');

xlim([min(gammaList), max(gammaList)]);
set(gca, 'XTick', gammaList);
set(gca, 'XTickLabel', arrayfun(@num2str, gammaList, 'UniformOutput', false));

if showReferenceLine
    legendText = sprintf('slope $=%g$', refSlope);

    legend({'complementarity defect', legendText}, ...
        'Interpreter', 'latex', ...
        'FontSize', 16, ...
        'Location', 'southwest');
else
    legend({'complementarity defect'}, ...
        'Interpreter', 'latex', ...
        'FontSize', 16, ...
        'Location', 'southwest');
end

% Save numerical values.
save(fullfile(outdir, sprintf('example3_complementarity_N%d.mat', Nfixed)), ...
    'gammaList', 'Nfixed', 'TuseList', 'comp_final');

% Save figure.
if saveFigure
    print(fig, fullfile(outdir, sprintf('example3_complementarity_N%d.eps', Nfixed)), '-depsc2');
    print(fig, fullfile(outdir, sprintf('example3_complementarity_N%d.png', Nfixed)), '-dpng', '-r300');
end
