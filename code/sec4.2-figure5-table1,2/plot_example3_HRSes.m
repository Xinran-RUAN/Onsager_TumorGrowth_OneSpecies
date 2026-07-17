% plot_example3_HSRes.m
% Compute and plot the final-time Hele--Shaw pressure residual
%
%   HSRes_{gamma,eta}(T)
%     = h * sum_{i in I_eta} | (delta_h^2 p_gamma^K)_i + 1 - p_gamma_i^K |,
%
% where
%
%   I_eta = { i : p_infty(x_i,T) >= eta }.
%
% This uses the same data files as plot_example3_comp.m.

clear; clc; close all;

outdir = 'example3';

gammaList = [10, 20, 40, 80, 160, 320];
Nfixed = 6400;

eta = 0.2;

showReferenceLine = true;
saveFigure = true;

HSRes_final = zeros(length(gammaList), 1);
numInteriorPts = zeros(length(gammaList), 1);
TuseList = zeros(length(gammaList), 1);
RexList = zeros(length(gammaList), 1);

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

    % Use load instead of matfile to avoid partial-loading warnings.
    S = load(matfileName);

    if ~isfield(S, 'sol_p')
        error('File %s does not contain sol_p.', files(idx).name);
    end

    % Mesh size.
    if isfield(S, 'h')
        h = S.h;
    elseif isfield(S, 'a') && isfield(S, 'b') && isfield(S, 'N')
        h = (S.b - S.a) / S.N;
    else
        h = 10 / Nfixed;  % fallback for the current [-5,5] test
    end

    % Grid points.
    if isfield(S, 'x')
        x = S.x(:);
    elseif isfield(S, 'a') && isfield(S, 'b') && isfield(S, 'N')
        a = S.a;
        N = S.N;
        x = a + ((1:N)' - 0.5) * h;
    else
        a = -5;
        x = a + ((1:Nfixed)' - 0.5) * h;
    end

    % Final-time pressure.
    sol_p = S.sol_p;
    if size(sol_p, 1) == length(x)
        p_final = sol_p(:, end);
    elseif size(sol_p, 2) == length(x)
        p_final = sol_p(end, :).';
    else
        error('Unexpected size of sol_p in %s.', files(idx).name);
    end
    p_final = p_final(:);

    % Final time.
    if isfield(S, 'time')
        Tuse = S.time(end);
    elseif isfield(S, 'T')
        Tuse = S.T;
    else
        error('File %s does not contain time or T.', files(idx).name);
    end
    TuseList(ig) = Tuse;

    % Exact Hele--Shaw radius.
    if isfield(S, 'R_exact')
        Rex = S.R_exact(end);
    elseif isfield(S, 'R0')
        Rex = asinh(exp(Tuse) * sinh(S.R0));
    else
        error('File %s does not contain R_exact or R0.', files(idx).name);
    end
    RexList(ig) = Rex;

    % Interior saturated set I_eta.
    % Also exclude the first and last grid points because delta_h^2 needs neighbors.
    p_inf_final = 1 - cosh(x) / cosh(Rex);
    p_inf_final = max(p_inf_final, 0);

    idxInterior = find(p_inf_final >= eta);
    idxInterior = idxInterior(idxInterior >= 2 & idxInterior <= length(x) - 1);

    if isempty(idxInterior)
        error('Interior set I_eta is empty for gamma = %d. Try smaller eta.', gamma);
    end

    % Second-order centered difference.
    d2p = (p_final(idxInterior + 1) ...
        - 2 * p_final(idxInterior) ...
        + p_final(idxInterior - 1)) / h^2;

    residual = d2p + 1 - p_final(idxInterior);

    HSRes_final(ig) = h * sum(abs(residual));
    numInteriorPts(ig) = length(idxInterior);

    fprintf('  T = %.8e, R(T) = %.8e, #I_eta = %d, HSRes = %.8e\n', ...
        Tuse, Rex, numInteriorPts(ig), HSRes_final(ig));
end

fprintf('\nSummary\n');
fprintf(' gamma        T              R(T)           #I_eta        HSRes_gamma_eta(T)\n');
for ig = 1:length(gammaList)
    fprintf('%6d   %12.6e   %12.6e   %8d      %12.6e\n', ...
        gammaList(ig), TuseList(ig), RexList(ig), ...
        numInteriorPts(ig), HSRes_final(ig));
end

% ========================================================================
% Reference line with slope -1, only as a visual guide.
% ========================================================================

if showReferenceLine
    refSlope = -1;
    iRef = 2;  % anchor at gamma = 20
    gamma_ref = gammaList(iRef);
    y_ref = 0.7 * HSRes_final(iRef);
    refLine = y_ref * (gammaList / gamma_ref).^refSlope;
end

% ========================================================================
% Plot
% ========================================================================

fig = figure;
set(fig, 'Color', 'w');
set(fig, 'Position', [100, 100, 700, 480]);

loglog(gammaList, HSRes_final, 'o-', ...
    'LineWidth', 2.2, ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', [0.00, 0.45, 0.74], ...
    'Color', [0.00, 0.45, 0.74]);
hold on;

if showReferenceLine
    loglog(gammaList, refLine, 'k--', 'LineWidth', 2.0);
end

grid on;
grid minor;
box on;

xlabel('$\gamma$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$\mathrm{HSRes}_{\gamma,\eta}(T)$', ...
    'Interpreter', 'latex', 'FontSize', 18);

title(sprintf('$T=1,\\ \\eta=%.2g$', eta), ...
    'Interpreter', 'latex', 'FontSize', 18);

set(gca, ...
    'FontSize', 18, ...
    'LineWidth', 1.0, ...
    'XScale', 'log', ...
    'YScale', 'log', ...
    'TickLabelInterpreter', 'latex');

xlim([min(gammaList), max(gammaList)]);
set(gca, 'XTick', gammaList);
set(gca, 'XTickLabel', arrayfun(@(g) ['$' num2str(g) '$'], ...
    gammaList, 'UniformOutput', false));

if showReferenceLine
    slopeLabel = sprintf('slope $=%g$', refSlope);
    legend({'pressure residual', slopeLabel}, ...
        'Interpreter', 'latex', ...
        'FontSize', 16, ...
        'Location', 'best');
else
    legend({'pressure residual'}, ...
        'Interpreter', 'latex', ...
        'FontSize', 16, ...
        'Location', 'best');
end

% Save numerical values.
save(fullfile(outdir, sprintf('example3_HSRes_eta%.2g_N%d.mat', eta, Nfixed)), ...
    'gammaList', 'Nfixed', 'eta', ...
    'TuseList', 'RexList', 'numInteriorPts', ...
    'HSRes_final');

% Save figure.
if saveFigure
    print(fig, fullfile(outdir, sprintf('example3_HSRes_eta%.2g_N%d.png', eta, Nfixed)), ...
        '-dpng', '-r300');
    print(fig, fullfile(outdir, sprintf('example3_HSRes_eta%.2g_N%d.eps', eta, Nfixed)), ...
        '-depsc2');
end
