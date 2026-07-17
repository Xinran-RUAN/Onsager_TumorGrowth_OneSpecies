clear; clc; close all;

outdir = 'example5';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

% ------------------------------------------------------------
% Domain and mesh
% ------------------------------------------------------------
a = -5;
b = 5;
Nx = 501;
Ny = 501;
h = (b - a) / (Nx - 1);

x = linspace(a, b, Nx);
y = linspace(a, b, Ny);
[X, Y] = meshgrid(x, y);

% ------------------------------------------------------------
% Parameters
% ------------------------------------------------------------
tau = 1e-3;
T = 2.0;
Nt = round(T / tau);
time = (0:Nt)' * tau;

gamma = 40;
pH = 1.0;
alpha = 1.0;

Gfun = @(p) alpha * (pH - p);
GpH = -alpha;

% ------------------------------------------------------------
% Initial data: annulus
% ------------------------------------------------------------
r = sqrt(X.^2 + Y.^2);

r_in = 0.8;
r_out = 2.0;
n0 = 0.8 * double(r >= r_in & r <= r_out);

n = n0;
p = n.^gamma;

initial_type = 'annulus';

opts.tol = 1e-10;
opts.maxit = 30;
opts.verbose = false;
opts.lineSearchMaxit = 30;

% ------------------------------------------------------------
% Snapshot saving
% ------------------------------------------------------------
% Save snapshots every 0.01 time units.
snapshotDt = 0.01;
snapshotTimes = 0:snapshotDt:T;

% Make sure final time is included.
if abs(snapshotTimes(end) - T) > 1e-14
    snapshotTimes = [snapshotTimes, T];
end

snapshotSteps = round(snapshotTimes / tau) + 1;
snapshotSteps = unique(snapshotSteps, 'stable');
snapshotSteps(snapshotSteps < 1) = 1;
snapshotSteps(snapshotSteps > Nt + 1) = Nt + 1;

snapshotTimes = (snapshotSteps - 1) * tau;
numSnap = length(snapshotSteps);

% Each snapshot is saved into a separate MAT-file.
snapshot_dir = fullfile(outdir, ...
    sprintf('snapshots_%s_gamma%d_N%d', initial_type, gamma, Nx));

if ~exist(snapshot_dir, 'dir')
    mkdir(snapshot_dir);
end

% Clean old snapshot files to avoid mixing different runs.
clearOldSnapshots = true;
if clearOldSnapshots
    oldFiles = dir(fullfile(snapshot_dir, 'snap_*.mat'));
    for jf = 1:length(oldFiles)
        delete(fullfile(snapshot_dir, oldFiles(jf).name));
    end
end

snap_files = cell(numSnap, 1);
snap_t = snapshotTimes(:);

% Save initial snapshot.
snap_id = 1;
tnow = 0.0;
step = 0;

snap_files{snap_id} = save_one_snapshot( ...
    snapshot_dir, snap_id, tnow, step, ...
    n, p, gamma, pH, alpha, initial_type);

fprintf('Saved snapshot %4d / %4d at t = %.4f\n', ...
    snap_id, numSnap, tnow);

snap_id = snap_id + 1;

% ------------------------------------------------------------
% Diagnostics
% ------------------------------------------------------------
min_n = zeros(Nt+1, 1);
max_n = zeros(Nt+1, 1);
max_p = zeros(Nt+1, 1);
energy = zeros(Nt+1, 1);
newton_iter = zeros(Nt, 1);
newton_res = zeros(Nt, 1);

min_n(1) = min(n(:));
max_n(1) = max(n(:));
max_p(1) = max(p(:));
energy(1) = h^2 * sum(n(:).^(gamma+1) / (gamma+1) - pH * n(:));

% ------------------------------------------------------------
% Plot during computation
% ------------------------------------------------------------
plotEvery = max(round(0.05 / tau), 1);

figure;
subplot(1,2,1);
imagesc(x, y, n);
axis equal tight;
set(gca, 'YDir', 'normal');
colorbar;
caxis([0, 1]);
title(sprintf('Density, t = %.3f', 0));

subplot(1,2,2);
imagesc(x, y, p);
axis equal tight;
set(gca, 'YDir', 'normal');
colorbar;
caxis([0, 1]);
title(sprintf('Pressure, t = %.3f', 0));
drawnow;

% ------------------------------------------------------------
% Time evolution
% ------------------------------------------------------------
for k = 1:Nt

    [n, p, info] = tumor_step_newton_2d(n, h, tau, gamma, pH, Gfun, GpH, opts);

    if ~info.converged
        fprintf('Step %d: Newton not fully converged, residual = %.4e\n', ...
                k, info.residual);
    end

    min_n(k+1) = min(n(:));
    max_n(k+1) = max(n(:));
    max_p(k+1) = max(p(:));
    energy(k+1) = h^2 * sum(n(:).^(gamma+1) / (gamma+1) - pH * n(:));

    newton_iter(k) = info.iter;
    newton_res(k) = info.residual;

    % --------------------------------------------------------
    % Save snapshots separately
    % --------------------------------------------------------
    while snap_id <= numSnap && k + 1 >= snapshotSteps(snap_id)

        tnow = k * tau;
        step = k;

        snap_files{snap_id} = save_one_snapshot( ...
            snapshot_dir, snap_id, tnow, step, ...
            n, p, gamma, pH, alpha, initial_type);

        fprintf('Saved snapshot %4d / %4d at t = %.4f\n', ...
            snap_id, numSnap, tnow);

        snap_id = snap_id + 1;
    end

    % --------------------------------------------------------
    % Online plotting
    % --------------------------------------------------------
    if mod(k, plotEvery) == 0 || k == Nt

        subplot(1,2,1);
        imagesc(x, y, n);
        axis equal tight;
        set(gca, 'YDir', 'normal');
        colorbar;
        caxis([0, 1]);
        title(sprintf('Density, t = %.3f', k*tau));

        subplot(1,2,2);
        imagesc(x, y, p);
        axis equal tight;
        set(gca, 'YDir', 'normal');
        colorbar;
        caxis([0, 1]);
        title(sprintf('Pressure, t = %.3f', k*tau));

        drawnow;
    end
end

% ------------------------------------------------------------
% Save summary data
% ------------------------------------------------------------
summaryfile = fullfile(outdir, ...
    sprintf('example5_2d_%s_gamma%d_N%d_summary.mat', initial_type, gamma, Nx));

save(summaryfile, ...
    'a', 'b', 'Nx', 'Ny', 'h', 'x', 'y', ...
    'tau', 'T', 'Nt', 'time', ...
    'gamma', 'pH', 'alpha', 'initial_type', ...
    'r_in', 'r_out', ...
    'snapshotDt', 'snapshotTimes', 'snapshotSteps', ...
    'snap_t', 'snapshot_dir', 'snap_files', ...
    'min_n', 'max_n', 'max_p', 'energy', ...
    'newton_iter', 'newton_res');

fprintf('\nExample 5 annulus finished.\n');
fprintf('Summary data saved to %s\n', summaryfile);
fprintf('Snapshots saved in %s\n', snapshot_dir);
fprintf('min n = %.6e, max n = %.6e, max p = %.6e\n', ...
        min_n(end), max_n(end), max_p(end));
fprintf('energy change = %.6e\n', energy(end) - energy(1));


% ========================================================================
% Local functions
% ========================================================================

function filename = save_one_snapshot(snapshot_dir, snap_id, tnow, step, ...
    n, p, gamma, pH, alpha, initial_type)
%SAVE_ONE_SNAPSHOT Save one snapshot into an individual MAT-file.

    ttag = time_tag(tnow);

    filename = fullfile(snapshot_dir, ...
        sprintf('snap_%04d_t_%s.mat', snap_id, ttag));

    save(filename, ...
        'n', 'p', 'tnow', 'step', ...
        'gamma', 'pH', 'alpha', 'initial_type');
end


function tag = time_tag(t)
%TIME_TAG Convert time value into a filename-friendly string.
%
% Example:
%   0       -> 0p000
%   0.01    -> 0p010
%   1.2     -> 1p200

    tag = sprintf('%.3f', t);
    tag = strrep(tag, '.', 'p');
    tag = strrep(tag, '-', 'm');
end
