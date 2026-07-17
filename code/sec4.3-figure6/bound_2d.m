% summarize_2d_structure_diagnostics.m
% Structure diagnostics for the two-dimensional examples.
%
% It computes
%   min n,
%   max n,
%   max p,
%   E_h(n^0),
%   E_h(n^K),
%   D_E^+ = max_k (E_h(n^{k+1}) - E_h(n^k))_+.
%
% The script first tries to read summary files. If no summary file is found,
% it computes the diagnostics from the stored snapshot files.

clear; clc;
format long

% ------------------------------------------------------------
% Basic settings
% ------------------------------------------------------------

gamma = 40;
pH = 1;

a = -5;
b = 5;
c = -5;
d = 5;

% ------------------------------------------------------------
% Cases
% ------------------------------------------------------------

cases = struct([]);

cases(1).name = 'annulus filling';
cases(1).summaryPatterns = { ...
    '*annulus*gamma40*N501*summary*.mat' ...
    };
cases(1).snapdirPatterns = { ...
    'snapshots*annulus*gamma40*N501*' ...
    };

cases(2).name = 'patch merging';
cases(2).summaryPatterns = { ...
    '*merging*gamma40*N501*summary*.mat', ...
    '*merge*gamma40*N501*summary*.mat', ...
    '*patch*gamma40*N501*summary*.mat', ...
    '*patches*gamma40*N501*summary*.mat', ...
    '*two*p*tch*gamma40*N501*summary*.mat', ...
    '*twodisks*gamma40*N501*summary*.mat', ...
    '*two*disk*gamma40*N501*summary*.mat' ...
    };
cases(2).snapdirPatterns = { ...
    'snapshots*merging*gamma40*N501*', ...
    'snapshots*merge*gamma40*N501*', ...
    'snapshots*patch*gamma40*N501*', ...
    'snapshots*patches*gamma40*N501*', ...
    'snapshots*two*p*tch*gamma40*N501*', ...
    'snapshots*twodisks*gamma40*N501*', ...
    'snapshots*two*disk*gamma40*N501*' ...
    };

% ------------------------------------------------------------
% Compute diagnostics
% ------------------------------------------------------------

numCases = length(cases);

minN = zeros(numCases, 1);
maxN = zeros(numCases, 1);
maxP = zeros(numCases, 1);
energyInitial = zeros(numCases, 1);
energyFinal = zeros(numCases, 1);
maxEnergyIncrease = zeros(numCases, 1);
sourceText = cell(numCases, 1);

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

fprintf('\nSearching data under: %s\n', searchRoot);

for ic = 1:numCases

    fprintf('\nCase: %s\n', cases(ic).name);

    [summaryFile, hasSummary] = find_file_recursive(searchRoot, cases(ic).summaryPatterns);

    if hasSummary
        fprintf('  Reading summary file: %s\n', summaryFile);
        S = load(summaryFile);

        [ok, minN(ic), maxN(ic), maxP(ic), ...
            energyInitial(ic), energyFinal(ic), maxEnergyIncrease(ic)] = ...
            diagnostics_from_summary(S);

        if ok
            sourceText{ic} = 'summary';
            continue;
        else
            fprintf('  Summary file found but necessary fields are incomplete.\n');
            fprintf('  Falling back to snapshots.\n');
        end
    end

    [snapdir, hasSnapdir] = find_dir_recursive(searchRoot, cases(ic).snapdirPatterns);

    if ~hasSnapdir
        error('No usable data found for case: %s.', cases(ic).name);
    end

    fprintf('  Reading snapshots from: %s\n', snapdir);
    fprintf('  Computing diagnostics from saved snapshots.\n');

    [minN(ic), maxN(ic), maxP(ic), ...
        energyInitial(ic), energyFinal(ic), maxEnergyIncrease(ic)] = ...
        diagnostics_from_snapshots(snapdir, gamma, pH, a, b, c, d);

    sourceText{ic} = 'saved snapshots';
end

% ------------------------------------------------------------
% Print summary
% ------------------------------------------------------------

fprintf('\nTwo-dimensional structure diagnostics\n');
fprintf('test                     min n          max n          max p          E_h(n^0)       E_h(n^K)       D_E^+          source\n');

for ic = 1:numCases
    fprintf('%-22s  %12.8e   %12.8e   %12.8e   %12.8e   %12.8e   %12s   %s\n', ...
        cases(ic).name, ...
        minN(ic), maxN(ic), maxP(ic), ...
        energyInitial(ic), energyFinal(ic), deplus_text(maxEnergyIncrease(ic)), ...
        sourceText{ic});
end

% ------------------------------------------------------------
% Write LaTeX table
% ------------------------------------------------------------

texfile = 'table_2d_structure_diagnostics.tex';
fid = fopen(texfile, 'w');

fprintf(fid, '\\begin{table}[t]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\caption{Structure diagnostics for the two-dimensional tests.}\n');
fprintf(fid, '\\label{tab:2d-structure}\n');
fprintf(fid, '\\begingroup\n');
fprintf(fid, '\\small\n');
fprintf(fid, '\\renewcommand{\\arraystretch}{1.08}\n');
fprintf(fid, '\\setlength{\\tabcolsep}{4pt}\n');
fprintf(fid, '\\begin{tabular}{lcccccc}\n');
fprintf(fid, '\\toprule\n');
fprintf(fid, 'test & $\\min n$ & $\\max n$ & $\\max p$ & $E_h(n^0)$ & $E_h(n^K)$ & $D_E^+$ \\\\\n');
fprintf(fid, '\\midrule\n');

for ic = 1:numCases
    fprintf(fid, '%s & %s & %s & %s & %s & %s & %s \\\\\n', ...
        cases(ic).name, ...
        sci_latex(minN(ic)), ...
        sci_latex(maxN(ic)), ...
        sci_latex(maxP(ic)), ...
        sci_latex(energyInitial(ic)), ...
        sci_latex(energyFinal(ic)), ...
        sci_latex_deplus(maxEnergyIncrease(ic)));
end

fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\endgroup\n');
fprintf(fid, '\\end{table}\n');

fclose(fid);

fprintf('\nSaved LaTeX table to %s\n', texfile);

% ========================================================================
% Local functions
% ========================================================================

function [fname, ok] = find_file_recursive(rootdir, patterns)
%FIND_FILE_RECURSIVE Find newest file matching any pattern under rootdir.

    fname = '';
    ok = false;

    allPaths = strsplit(genpath(rootdir), pathsep);

    bestDate = -inf;
    bestFile = '';

    for ip = 1:length(allPaths)

        thisDir = allPaths{ip};
        if isempty(thisDir)
            continue;
        end

        for jp = 1:length(patterns)

            files = dir(fullfile(thisDir, patterns{jp}));

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

    if ~isempty(bestFile)
        fname = bestFile;
        ok = true;
    end
end


function [dname, ok] = find_dir_recursive(rootdir, patterns)
%FIND_DIR_RECURSIVE Find newest directory matching any pattern under rootdir.

    dname = '';
    ok = false;

    allPaths = strsplit(genpath(rootdir), pathsep);

    bestDate = -inf;
    bestDir = '';

    for ip = 1:length(allPaths)

        thisDir = allPaths{ip};
        if isempty(thisDir)
            continue;
        end

        for jp = 1:length(patterns)

            dirs = dir(fullfile(thisDir, patterns{jp}));

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

    if ~isempty(bestDir)
        dname = bestDir;
        ok = true;
    end
end


function [ok, minN, maxN, maxP, energyInitial, energyFinal, maxEnergyIncrease] = ...
    diagnostics_from_summary(S)

    ok = false;

    minN = NaN;
    maxN = NaN;
    maxP = NaN;
    energyInitial = NaN;
    energyFinal = NaN;
    maxEnergyIncrease = NaN;

    % Common possible names.
    minNArray = get_field_if_exists(S, {'min_n', 'minN', 'n_min'});
    maxNArray = get_field_if_exists(S, {'max_n', 'maxN', 'n_max'});
    maxPArray = get_field_if_exists(S, {'max_p', 'maxP', 'p_max'});
    energyArray = get_field_if_exists(S, {'energy', 'Eh', 'E', 'energy_hist'});

    if isempty(minNArray) || isempty(maxNArray) || isempty(maxPArray) || isempty(energyArray)
        return;
    end

    minN = min(minNArray(:));
    maxN = max(maxNArray(:));
    maxP = max(maxPArray(:));

    energyArray = energyArray(:);
    energyInitial = energyArray(1);
    energyFinal = energyArray(end);

    if length(energyArray) >= 2
        maxEnergyIncrease = max(diff(energyArray));
        maxEnergyIncrease = max(maxEnergyIncrease, 0);
    else
        maxEnergyIncrease = NaN;
    end

    ok = true;
end


function value = get_field_if_exists(S, names)

    value = [];

    for k = 1:length(names)
        if isfield(S, names{k})
            value = S.(names{k});
            return;
        end
    end
end


function [minN, maxN, maxP, energyInitial, energyFinal, maxEnergyIncrease] = ...
    diagnostics_from_snapshots(snapdir, gamma, pH, a, b, c, d)

    files = dir(fullfile(snapdir, 'snap_*.mat'));

    if isempty(files)
        error('No snapshot files found in %s.', snapdir);
    end

    [~, idx] = sort({files.name});
    files = files(idx);

    numFiles = length(files);
    snap_t = zeros(numFiles, 1);

    for k = 1:numFiles
        snap_t(k) = parse_time_from_filename(files(k).name);
    end

    [snap_t, idx] = sort(snap_t);
    files = files(idx);

    minN = inf;
    maxN = -inf;
    maxP = -inf;

    energy = zeros(numFiles, 1);

    for k = 1:numFiles

        S = load(fullfile(snapdir, files(k).name));
        [n, p, ~] = get_snapshot_vars(S);

        if isempty(p)
            p = n.^gamma;
        end

        [Ny, Nx] = size(n);
        hx = (b - a) / Nx;
        hy = (d - c) / Ny;

        minN = min(minN, min(n(:)));
        maxN = max(maxN, max(n(:)));
        maxP = max(maxP, max(p(:)));

        energy(k) = hx * hy * sum( ...
            1/(gamma+1) * n(:).^(gamma+1) - pH * n(:) );
    end

    energyInitial = energy(1);
    energyFinal = energy(end);

    if length(energy) >= 2
        maxEnergyIncrease = max(diff(energy));
        maxEnergyIncrease = max(maxEnergyIncrease, 0);
    else
        maxEnergyIncrease = NaN;
    end
end


function t = parse_time_from_filename(fname)

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


function [n, p, tnow] = get_snapshot_vars(S)

    if isfield(S, 'n')
        n = S.n;
    elseif isfield(S, 'n_snap')
        n = S.n_snap;
    else
        error('Snapshot does not contain n or n_snap.');
    end

    if isfield(S, 'p')
        p = S.p;
    elseif isfield(S, 'p_snap')
        p = S.p_snap;
    else
        p = [];
    end

    if isfield(S, 'tnow')
        tnow = S.tnow;
    elseif isfield(S, 't_snap')
        tnow = S.t_snap;
    else
        tnow = NaN;
    end
end


function s = sci_latex(x)

    if isnan(x)
        s = '--';
        return;
    end

    if abs(x) < 1e-14
        s = '$0$';
        return;
    end

    exponent = floor(log10(abs(x)));
    mantissa = x / 10^exponent;

    s = sprintf('$%.3f{\\rm E}{%+d}$', mantissa, exponent);
    s = strrep(s, '{+','{');
end


function s = sci_latex_deplus(x)

    if isnan(x)
        s = '--';
        return;
    end

    if x < 1e-14
        s = '$\le 10^{-14}$';
        return;
    end

    s = sci_latex(x);
end


function s = deplus_text(x)

    if isnan(x)
        s = 'NaN';
        return;
    end

    if x < 1e-14
        s = '<= 1e-14';
        return;
    end

    s = sprintf('%.8e', x);
end
