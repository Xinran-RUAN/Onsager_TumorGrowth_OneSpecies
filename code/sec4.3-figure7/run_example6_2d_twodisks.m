% run_example6_2d_twodisks.m
% Two-dimensional tumor-growth example: merging of two disks.
%
% Snapshots are saved every 0.01 time units as separate MAT files.

clear; clc; close all;

outdir = 'example6';
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
T = 4.0;
Nt = round(T / tau);
time = (0:Nt)' * tau;

gamma = 40;
pH = 1.0;
alpha = 1.0;

Gfun = @(p) alpha * (pH - p);
GpH = -alpha;

% ------------------------------------------------------------
% Initial data: two disks
% ------------------------------------------------------------
r0 = 0.8;
c0 = 1.4;

disk1 = (X + c0).^2 + Y.^2 <= r0^2;
disk2 = (X - c0).^2 + Y.^2 <= r0^2;

n0 = 0.8 * double(disk1 | disk2);
n = n0;
p = n.^gamma;

initial_type = 'twodisks';

opts.tol = 1e-9;
opts.maxit = 30;
opts.verbose = false;
opts.lineSearchMaxit = 30;

% ------------------------------------------------------------
% Save snapshots every 0.01 as separate files
% ------------------------------------------------------------
snapshotDt = 0.01;
snapshotTimes = 0:snapshotDt:T;

if abs(snapshotTimes(end) - T) > 1e-14
    snapshotTimes = [snapshotTimes, T];
end

snapshotSteps = round(snapshotTimes / tau) + 1;
snapshotSteps = unique(snapshotSteps, 'stable');
snapshotSteps(snapshotSteps < 1) = 1;
snapshotSteps(snapshotSteps > Nt + 1) = Nt + 1;

snapshotTimes = (snapshotSteps - 1) * tau;
numSnap = length(snapshotSteps);
snap_t = snapshotTimes(:);

case_name = sprintf('%s_gamma%d_N%d', initial_type, gamma, Nx);
snapshot_dir = fullfile(outdir, ['snapshots_', case_name]);

if ~exist(snapshot_dir, 'dir')
    mkdir(snapshot_dir);
end

% Remove old snapshot files to avoid mixing different runs.
clearOldSnapshots = true;
if clearOldSnapshots
    oldFiles = dir(fullfile(snapshot_dir, 'snap_*.mat'));
    for jf = 1:length(oldFiles)
        delete(fullfile(snapshot_dir, oldFiles(jf).name));
    end
end

snap_files = cell(numSnap, 1);
snap_names = cell(numSnap, 1);

% Save initial snapshot.
snap_id = 1;
tnow = 0.0;
step = 0;

[snap_files{snap_id}, snap_names{snap_id}] = save_one_snapshot( ...
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

figRun = figure;
set(figRun, 'Color', 'w');

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

        step = k;
        tnow = k * tau;

        [snap_files{snap_id}, snap_names{snap_id}] = save_one_snapshot( ...
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
datafile = fullfile(outdir, ...
    sprintf('example6_2d_%s_gamma%d_N%d_summary.mat', initial_type, gamma, Nx));

snapdir = snapshot_dir;

save(datafile, ...
    'a', 'b', 'Nx', 'Ny', 'h', 'x', 'y', 'X', 'Y', ...
    'tau', 'T', 'Nt', 'time', ...
    'gamma', 'pH', 'alpha', 'initial_type', ...
    'r0', 'c0', ...
    'snapshotDt', 'snapshotTimes', 'snapshotSteps', ...
    'snap_t', 'snapshot_dir', 'snapdir', 'snap_files', 'snap_names', ...
    'min_n', 'max_n', 'max_p', 'energy', ...
    'newton_iter', 'newton_res');

fprintf('\nExample 6 two disks finished.\n');
fprintf('Summary data saved to %s\n', datafile);
fprintf('Snapshots saved in %s\n', snapshot_dir);
fprintf('min n = %.6e, max n = %.6e, max p = %.6e\n', ...
        min_n(end), max_n(end), max_p(end));
fprintf('energy change = %.6e\n', energy(end) - energy(1));

% If you want to plot immediately after computation, uncomment the line below.
% plot_example6_twodisks;


% ========================================================================
% Local functions
% ========================================================================

function [fullfile_name, short_name] = save_one_snapshot( ...
    snapshot_dir, snap_id, tnow, step, ...
    n, p, gamma, pH, alpha, initial_type)
%SAVE_ONE_SNAPSHOT Save one snapshot into an individual MAT-file.

    short_name = make_snapshot_name(snap_id, tnow);
    fullfile_name = fullfile(snapshot_dir, short_name);

    save(fullfile_name, ...
        'n', 'p', 'tnow', 'step', ...
        'gamma', 'pH', 'alpha', 'initial_type');
end


function fname = make_snapshot_name(snap_id, tnow)
%MAKE_SNAPSHOT_NAME Generate a filename-friendly snapshot name.

    tstr = sprintf('%.4f', tnow);
    tstr = strrep(tstr, '.', 'p');
    tstr = strrep(tstr, '-', 'm');

    fname = sprintf('snap_%04d_t_%s.mat', snap_id, tstr);
end