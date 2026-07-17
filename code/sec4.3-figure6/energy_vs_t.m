% plot_2d_energy_evolution.m
% Plot modified energy evolution for the two-dimensional tests.
%
% x-axis: time
% y-axis: E_h(n^k)
%
% This script reads summary files, for example:
%   example5_2d_annulus_gamma40_N501_summary.mat
%   example6_2d_twodisks_gamma40_N501_summary.mat

clear; clc; close all;

% ------------------------------------------------------------
% Settings
% ------------------------------------------------------------

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end
[~, scriptFolder] = fileparts(scriptDir);
if length(scriptFolder) >= 3 && strcmp(scriptFolder(end-2:end), '_2D')
    searchRoot = fileparts(scriptDir);
else
    searchRoot = scriptDir;
end

saveFigure = true;
outfilePrefix = 'Figure_2d_energy_evolution';

% ------------------------------------------------------------
% Locate data files
% ------------------------------------------------------------

annulusPatterns = { ...
    '*annulus*gamma40*N501*summary*.mat' ...
    };
annulusSnapdirPatterns = { ...
    'snapshots*annulus*gamma40*N501*' ...
    };

mergingPatterns = { ...
    '*twodisks*gamma40*N501*summary*.mat', ...
    '*two*disk*gamma40*N501*summary*.mat', ...
    '*merging*gamma40*N501*summary*.mat', ...
    '*merge*gamma40*N501*summary*.mat', ...
    '*patch*gamma40*N501*summary*.mat', ...
    '*patches*gamma40*N501*summary*.mat' ...
    };
mergingSnapdirPatterns = { ...
    'snapshots*twodisks*gamma40*N501*', ...
    'snapshots*two*disk*gamma40*N501*', ...
    'snapshots*merging*gamma40*N501*', ...
    'snapshots*merge*gamma40*N501*', ...
    'snapshots*patch*gamma40*N501*', ...
    'snapshots*patches*gamma40*N501*' ...
    };

annulusFile = find_file_recursive(searchRoot, annulusPatterns);
mergingFile = find_file_recursive(searchRoot, mergingPatterns);
annulusSnapdir = find_dir_recursive(searchRoot, annulusSnapdirPatterns);
mergingSnapdir = find_dir_recursive(searchRoot, mergingSnapdirPatterns);

fprintf('\nSearching summary files under: %s\n', searchRoot);

if isempty(annulusFile)
    error('Annulus summary file not found.');
end

if isempty(mergingFile)
    error('Merging/two-disks summary file not found.');
end

fprintf('Annulus file: %s\n', annulusFile);
fprintf('Merging file: %s\n', mergingFile);
fprintf('Annulus snapshots: %s\n', annulusSnapdir);
fprintf('Merging snapshots: %s\n', mergingSnapdir);

% ------------------------------------------------------------
% Read data
% ------------------------------------------------------------

[tAnn, EAnn] = read_time_energy_with_snapshot_fallback( ...
    annulusFile, annulusSnapdir, 'annulus filling');
[tMer, EMer] = read_time_energy_with_snapshot_fallback( ...
    mergingFile, mergingSnapdir, 'patch merging');

fprintf('Annulus time interval: [%.4g, %.4g]\n', min(tAnn), max(tAnn));
fprintf('Merging time interval: [%.4g, %.4g]\n', min(tMer), max(tMer));

% ------------------------------------------------------------
% Plot
% ------------------------------------------------------------

fig = figure;
set(fig, 'Color', 'w');
set(fig, 'Position', [100, 100, 760, 500]);

hold on;

plot(tAnn, EAnn, '-', ...
    'LineWidth', 2.4);

plot(tMer, EMer, '--', ...
    'LineWidth', 2.4);


box on;
grid on;

xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$E_h(n^k)$', 'Interpreter', 'latex', 'FontSize', 20);
legend({'annulus filling', 'patch merging'}, ...
    'Interpreter', 'latex', ...
    'FontSize', 16, ...
    'Location', 'best');

set(gca, ...
    'FontSize', 18, ...
    'LineWidth', 1.0, ...
    'TickLabelInterpreter', 'latex');

xlim([0, max([tAnn(:); tMer(:)])]);

% ------------------------------------------------------------
% Check monotonicity
% ------------------------------------------------------------

dEAnn = diff(EAnn);
dEMer = diff(EMer);

fprintf('\nEnergy monotonicity check\n');
fprintf('annulus filling: max positive increment = %.8e\n', max(max(dEAnn), 0));
fprintf('patch merging:   max positive increment = %.8e\n', max(max(dEMer), 0));

% ------------------------------------------------------------
% Save figure
% ------------------------------------------------------------

if saveFigure
    print(fig, [outfilePrefix, '.png'], '-dpng', '-r300');
    print(fig, [outfilePrefix, '.eps'], '-depsc2');
    savefig(fig, [outfilePrefix, '.fig']);

    fprintf('\nSaved figure to %s.png\n', outfilePrefix);
end

% ========================================================================
% Local functions
% ========================================================================

function fname = find_file_recursive(rootdir, patterns)
%FIND_FILE_RECURSIVE Find newest file matching any pattern under rootdir.

    fname = '';

    allDirs = strsplit(genpath(rootdir), pathsep);

    bestDate = -inf;
    bestFile = '';

    for id = 1:length(allDirs)

        thisDir = allDirs{id};
        if isempty(thisDir)
            continue;
        end

        for ip = 1:length(patterns)

            files = dir(fullfile(thisDir, patterns{ip}));

            for k = 1:length(files)

                if files(k).isdir
                    continue;
                end

                if files(k).datenum > bestDate
                    bestDate = files(k).datenum;
                    bestFile = fullfile(thisDir, files(k).name);
                end
            end
        end
    end

    fname = bestFile;
end


function dname = find_dir_recursive(rootdir, patterns)
%FIND_DIR_RECURSIVE Find newest directory matching any pattern under rootdir.

    dname = '';

    allDirs = strsplit(genpath(rootdir), pathsep);

    bestDate = -inf;
    bestDir = '';

    for id = 1:length(allDirs)

        thisDir = allDirs{id};
        if isempty(thisDir)
            continue;
        end

        for ip = 1:length(patterns)

            dirs = dir(fullfile(thisDir, patterns{ip}));

            for k = 1:length(dirs)

                if ~dirs(k).isdir
                    continue;
                end

                name = dirs(k).name;
                if strcmp(name, '.') || strcmp(name, '..')
                    continue;
                end

                if dirs(k).datenum > bestDate
                    bestDate = dirs(k).datenum;
                    bestDir = fullfile(thisDir, name);
                end
            end
        end
    end

    dname = bestDir;
end


function [t, E] = read_time_energy_with_snapshot_fallback(filename, snapdir, label)
%READ_TIME_ENERGY_WITH_SNAPSHOT_FALLBACK Prefer summary, repair stale summary.

    [t, E] = read_time_energy(filename);

    if isempty(snapdir)
        return;
    end

    snapTimes = read_snapshot_times(snapdir);
    if isempty(snapTimes)
        return;
    end

    if max(snapTimes) > max(t) + 1e-12
        fprintf(['  %s summary ends at t = %.4g, but saved snapshots ', ...
            'extend to t = %.4g.\n'], label, max(t), max(snapTimes));
        fprintf('  Recomputing %s energy from saved snapshots.\n', label);
        [t, E] = read_snapshot_energy(snapdir, filename);
    end
end


function [t, E] = read_time_energy(filename)
%READ_TIME_ENERGY Read time and energy arrays from a summary file.

    S = load(filename);

    % Read energy.
    if isfield(S, 'energy')
        E = S.energy(:);
    elseif isfield(S, 'Eh')
        E = S.Eh(:);
    elseif isfield(S, 'E')
        E = S.E(:);
    elseif isfield(S, 'energy_hist')
        E = S.energy_hist(:);
    else
        error('File %s does not contain an energy array.', filename);
    end

    % Read time.
    if isfield(S, 'time')
        t = S.time(:);
    elseif isfield(S, 't')
        t = S.t(:);
    elseif isfield(S, 't_hist')
        t = S.t_hist(:);
    elseif isfield(S, 'tau')
        t = (0:length(E)-1)' * S.tau;
    elseif isfield(S, 'T')
        t = linspace(0, S.T, length(E)).';
    else
        warning('No time array found in %s. Using index as time.', filename);
        t = (0:length(E)-1).';
    end

    % Align lengths if needed.
    if length(t) == length(E) + 1
        t = t(1:end-1);
    elseif length(E) == length(t) + 1
        E = E(1:end-1);
    end

    if length(t) ~= length(E)
        error('Time and energy lengths do not match in %s.', filename);
    end
end


function snapTimes = read_snapshot_times(snapdir)
%READ_SNAPSHOT_TIMES Read physical times encoded in snapshot filenames.

    files = dir(fullfile(snapdir, 'snap_*.mat'));
    snapTimes = zeros(length(files), 1);

    for k = 1:length(files)
        snapTimes(k) = parse_time_from_filename(files(k).name);
    end

    snapTimes = snapTimes(isfinite(snapTimes));
end


function [t, E] = read_snapshot_energy(snapdir, summaryFile)
%READ_SNAPSHOT_ENERGY Recompute modified energy from saved snapshots.

    Ssummary = load(summaryFile);

    gamma = get_scalar_field(Ssummary, 'gamma', 40);
    pH = get_scalar_field(Ssummary, 'pH', 1);

    files = dir(fullfile(snapdir, 'snap_*.mat'));
    if isempty(files)
        error('No snapshot files found in %s.', snapdir);
    end

    snapTimes = zeros(length(files), 1);
    for k = 1:length(files)
        snapTimes(k) = parse_time_from_filename(files(k).name);
    end

    [t, idx] = sort(snapTimes);
    files = files(idx);

    E = zeros(length(files), 1);

    for k = 1:length(files)
        S = load(fullfile(snapdir, files(k).name));
        [n, tnow] = get_snapshot_vars(S);

        if isfinite(tnow)
            t(k) = tnow;
        end

        if isfield(Ssummary, 'h')
            hx = Ssummary.h;
            hy = Ssummary.h;
        else
            [Ny, Nx] = size(n);
            a = get_scalar_field(Ssummary, 'a', -5);
            b = get_scalar_field(Ssummary, 'b', 5);
            c = get_scalar_field(Ssummary, 'c', -5);
            d = get_scalar_field(Ssummary, 'd', 5);
            hx = (b - a) / Nx;
            hy = (d - c) / Ny;
        end

        E(k) = hx * hy * sum( ...
            1/(gamma+1) * n(:).^(gamma+1) - pH * n(:));
    end
end


function t = parse_time_from_filename(fname)
%PARSE_TIME_FROM_FILENAME Parse time from a snapshot filename.

    token = regexp(fname, '_t_([0-9mp]+)\.mat', 'tokens', 'once');

    if isempty(token)
        t = NaN;
        return;
    end

    s = token{1};
    s = strrep(s, 'p', '.');
    s = strrep(s, 'm', '-');

    t = str2double(s);
end


function [n, tnow] = get_snapshot_vars(S)
%GET_SNAPSHOT_VARS Read density and time from a snapshot structure.

    if isfield(S, 'n')
        n = S.n;
    elseif isfield(S, 'n_snap')
        n = S.n_snap;
    else
        error('Snapshot does not contain n or n_snap.');
    end

    if isfield(S, 'tnow')
        tnow = S.tnow;
    elseif isfield(S, 't_snap')
        tnow = S.t_snap;
    else
        tnow = NaN;
    end
end


function value = get_scalar_field(S, name, defaultValue)
%GET_SCALAR_FIELD Read scalar field with a default.

    if isfield(S, name)
        value = S.(name);
        value = value(1);
    else
        value = defaultValue;
    end
end
