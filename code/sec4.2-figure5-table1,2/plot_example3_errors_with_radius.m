% plot_example3_errors_with_radius.m
% Plot final-time density, pressure, and radius errors with respect to
% the explicit Hele--Shaw limit for a fixed spatial mesh size N.
%
% The radius error is computed from the mass-based numerical radius
%     R_h(T) = 0.5 * h * sum_i n_i(T),
% which is natural for the symmetric one-dimensional patch limit.

clear; clc; close all;

outdir = 'example3';

gammaList = [10, 20, 40, 80, 160, 320];
Nfixed = 6400;

err_n_final = zeros(length(gammaList), 1);
err_p_final = zeros(length(gammaList), 1);
err_R_final = zeros(length(gammaList), 1);
R_num_final = zeros(length(gammaList), 1);
R_ex_final  = zeros(length(gammaList), 1);
T_final     = zeros(length(gammaList), 1);

for ig = 1:length(gammaList)

    gamma = gammaList(ig);

    % Only read files with the fixed N.
    filePattern = sprintf('example3_gamma%d_N%d*.mat', gamma, Nfixed);
    files = dir(fullfile(outdir, filePattern));

    if isempty(files)
        error('No data file found for gamma = %d and N = %d.', gamma, Nfixed);
    end

    % If multiple files exist, use the newest one.
    [~, idx] = max([files.datenum]);
    matfile = fullfile(outdir, files(idx).name);

    fprintf('Reading gamma = %d, N = %d from %s\n', ...
        gamma, Nfixed, files(idx).name);

    S = load(matfile);

    % Final-time L1 errors.
    if ~isfield(S, 'err_n_L1') || ~isfield(S, 'err_p_L1')
        error('File %s does not contain err_n_L1 or err_p_L1.', files(idx).name);
    end
    err_n_final(ig) = S.err_n_L1(end);
    err_p_final(ig) = S.err_p_L1(end);

    % Final time.
    if isfield(S, 'time')
        Tuse = S.time(end);
    elseif isfield(S, 'T')
        Tuse = S.T;
    else
        error('File %s does not contain time or T.', files(idx).name);
    end
    T_final(ig) = Tuse;

    % Mass-based numerical radius:
    % for the limiting symmetric patch, int n dx = 2 R(t).
    if ~isfield(S, 'sol_n')
        error('File %s does not contain sol_n, which is needed for the radius error.', files(idx).name);
    end
    if isfield(S, 'h')
        h = S.h;
    elseif isfield(S, 'a') && isfield(S, 'b') && isfield(S, 'N')
        h = (S.b - S.a) / S.N;
    else
        error('File %s does not contain h or domain information.', files(idx).name);
    end

    n_final = S.sol_n(:, end);
    R_num_final(ig) = 0.5 * h * sum(n_final);

    % Exact Hele--Shaw radius.
    if isfield(S, 'R_exact')
        R_ex_final(ig) = S.R_exact(end);
    elseif isfield(S, 'R0')
        R_ex_final(ig) = asinh(exp(Tuse) * sinh(S.R0));
    else
        error('File %s does not contain R_exact or R0.', files(idx).name);
    end

    err_R_final(ig) = abs(R_num_final(ig) - R_ex_final(ig));

    fprintf('  density error  = %.8e\n', err_n_final(ig));
    fprintf('  pressure error = %.8e\n', err_p_final(ig));
    fprintf('  radius error   = %.8e  (R_h = %.8e, R = %.8e)\n', ...
        err_R_final(ig), R_num_final(ig), R_ex_final(ig));
end

fprintf('\nSummary\n');
fprintf(' gamma        err_n_final        err_p_final        err_R_final\n');
for ig = 1:length(gammaList)
    fprintf('%6d    %12.6e    %12.6e    %12.6e\n', ...
        gammaList(ig), err_n_final(ig), err_p_final(ig), err_R_final(ig));
end

% ========================================================================
% Build a reference line with slope -1.
% This line is only a visual guide.
% ========================================================================

iRef = 2;   % use gammaList(2)=20 as the anchor point
gamma_ref = gammaList(iRef);
refSlope = -1;

% Put the reference line below the main curves.
y_ref = 0.2 * min([err_n_final(iRef), err_p_final(iRef), err_R_final(iRef)]);
refLine = y_ref * (gammaList / gamma_ref).^refSlope;

% ========================================================================
% Plot
% ========================================================================

fig = figure;
set(fig, 'Color', 'w');
set(fig, 'Position', [100, 100, 700, 480]);

loglog(gammaList, err_n_final, 'o-', ...
    'LineWidth', 2.2, ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', [0.00, 0.45, 0.74]);
hold on;

loglog(gammaList, err_p_final, 's-', ...
    'LineWidth', 2.2, ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', [0.85, 0.33, 0.10]);

loglog(gammaList, err_R_final, '^-', ...
    'LineWidth', 2.2, ...
    'MarkerSize', 8, ...
    'MarkerFaceColor', [0.47, 0.67, 0.19]);

loglog(gammaList, refLine, 'k--', ...
    'LineWidth', 1.8);

xlabel('$\gamma$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('error', 'Interpreter', 'latex', 'FontSize', 20);

title('$T=1$', 'Interpreter', 'latex', 'FontSize', 20);

slopeLabel = sprintf('slope $=%g$', refSlope);
legend({'density error', 'pressure error', 'radius error', slopeLabel}, ...
    'Interpreter', 'latex', ...
    'FontSize', 18, ...
    'Location', 'best');

set(gca, ...
    'FontSize', 18, ...
    'LineWidth', 1.0, ...
    'TickLabelInterpreter', 'latex');

box on;
grid on;

xlim([min(gammaList), max(gammaList)]);
set(gca, 'XTick', gammaList);
set(gca, 'XTickLabel', arrayfun(@(g) ['$' num2str(g) '$'], ...
    gammaList, 'UniformOutput', false));

% Make y-axis range cleaner.
allErr = [err_n_final(:); err_p_final(:); err_R_final(:); refLine(:)];
allErr = allErr(isfinite(allErr) & allErr > 0);

if ~isempty(allErr)
    ymin = min(allErr);
    ymax = max(allErr);
    if ymin < ymax
        ylim([10^floor(log10(ymin)), 10^ceil(log10(ymax))]);
    end
end

outfile = sprintf('Figure10_example3_errors_with_radius_N%d.png', Nfixed);
print(fig, fullfile(outdir, outfile), '-dpng', '-r300');

fprintf('\nSaved figure to %s\n', fullfile(outdir, outfile));

% Save the plotted data.
datafile = sprintf('example3_errors_with_radius_N%d.mat', Nfixed);
save(fullfile(outdir, datafile), ...
    'gammaList', 'Nfixed', ...
    'err_n_final', 'err_p_final', 'err_R_final', ...
    'R_num_final', 'R_ex_final', 'T_final', ...
    'gamma_ref', 'y_ref', 'refSlope', 'refLine');

fprintf('Saved error data to %s\n', fullfile(outdir, datafile));
