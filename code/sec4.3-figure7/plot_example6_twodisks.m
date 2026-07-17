% plot_example6_twodisks.m
% Plotting script for the annulus example in 2D.
%
% Compatible with:
%
%   1. New separated snapshot format:
%        example5_2d_annulus_gamma*_N*_summary.mat
%        snapshots_annulus_gamma*_N*/snap_*.mat
%
%   2. Separated snapshots without summary file:
%        snapshots_annulus_gamma*_N*/snap_*.mat
%      In this case the grid is inferred from the snapshot size and the
%      default domain [-5,5]^2.
%
%   3. Old all-in-one format:
%        example5_2d_annulus_gamma*_N*.mat
%      containing snap_n, snap_p, snap_t.
%
% It generates:
%   Figure 0: 2 x 3 top-view density snapshots
%   Figure 1: 2 x 3 density cross-sections along y = 0
%   Figure 2: 2 x 3 pressure cross-sections along y = 0
%
% The same six times are used in all figures.
% plot_example5_annulus.m
% Simple plotting script for the annulus example.
%
% It reads separated snapshot files directly from:
%   snapshots_annulus_gamma40_N501
%
% It generates:
%   Figure 1: top-view density snapshots
%   Figure 2: density cross-sections along y = 0

clear; clc; close all;

% ------------------------------------------------------------
% Basic settings
% ------------------------------------------------------------

% If this script is in the parent folder, use example5/snapshots...
% If this script is inside example5, use snapshots...
snapdir = 'snapshots_annulus_gamma40_N501';

if ~exist(snapdir, 'dir')
    snapdir = fullfile('example6', 'snapshots_twodisks_gamma40_N501');
end

if ~exist(snapdir, 'dir')
    error('Snapshot folder not found.');
end

% Output folder
[outdir, ~, ~] = fileparts(snapdir);
if isempty(outdir)
    outdir = '.';
end

% Domain
a = -5;
b = 5;

% Six times for both figures
desiredTimes = [0, .80, 1.46, 1.50, 2., 3];

% Plot windows
xwin2D = [-4, 4];
ywin2D = [-4, 4];

xwinProfile = [-4, 4];
nylim = [-0.05, 1.05];

% Free-boundary level
nLevel = 0.5;

% Ticks
xtickVals2D = -4:2:4;
ytickVals2D = -4:2:4;

xtickValsProfile = -4:1:4;
ytickValsProfile = 0:0.2:1.2;

% Plot style parameters
topViewFigurePosition = [50, 50, 1200, 700];
sectionFigurePosition = [50, 50, 1500, 700];

axisLabelFontSize = 22;
titleFontSize = 22;
tickFontSize = 20;
colorbarFontSize = 22;
colorbarLabelFontSize = 22;
superTitleFontSize = 22;

axisLineWidth = 1.2;
contourLineWidth = 1.5;
profileLineWidth = 2.5;
thresholdLineWidth = 1.2;
zeroLineWidth = 0.8;

topViewLeftShift = 0.035;
topViewRowGapShift = 0.015;
sectionColumnShift = [-0.035, 0.000, 0.035];
sectionLeftShift = -0.03;
sectionWidthScale = 1.1;

densityColorLimits = [0, 1];
colorbarPosition = [0.9, 0.15, 0.015, 0.70];
colorbarTicks = 0:0.2:1;
saveResolution = 300;

xLabelText = '$x$';
yLabelText = '$y$';
densitySectionLabelText = '$n(x,0,t)$';
colorbarLabelText = '$n$';
timeTitleFormat = '$t=%.2f$';

% ------------------------------------------------------------
% Read snapshot list
% ------------------------------------------------------------

files = dir(fullfile(snapdir, 'snap_*.mat'));

if isempty(files)
    error('No snapshot files found in %s.', snapdir);
end

% Sort by filename
[~, idx] = sort({files.name});
files = files(idx);

numFiles = length(files);
snap_t = zeros(numFiles, 1);

for k = 1:numFiles
    snap_t(k) = parse_time_from_filename(files(k).name);
end

% Sort by physical time
[snap_t, idx] = sort(snap_t);
files = files(idx);

fprintf('Snapshot folder: %s\n', snapdir);
fprintf('Number of snapshots = %d\n', numFiles);
fprintf('Time range = [%.4f, %.4f]\n', min(snap_t), max(snap_t));

if max(desiredTimes) > max(snap_t) + 1e-12
    error('Requested final time %.4f exceeds largest stored time %.4f.', ...
        max(desiredTimes), max(snap_t));
end

% ------------------------------------------------------------
% Load selected snapshots
% ------------------------------------------------------------

plotIds = zeros(size(desiredTimes));

for k = 1:length(desiredTimes)
    [~, plotIds(k)] = min(abs(snap_t - desiredTimes(k)));
end

nPanels = length(plotIds);

% Read first selected snapshot to determine grid size
S0 = load(fullfile(snapdir, files(plotIds(1)).name));
[n0, ~, ~] = get_snapshot_vars(S0);

[Ny, Nx] = size(n0);
x = linspace(a, b, Nx);
y = linspace(a, b, Ny);

[~, iy0] = min(abs(y));

nCell = cell(nPanels, 1);
tCell = zeros(nPanels, 1);

fprintf('\nSelected snapshots:\n');

for k = 1:nPanels
    id = plotIds(k);
    S = load(fullfile(snapdir, files(id).name));

    [nNow, ~, tNow] = get_snapshot_vars(S);

    if isnan(tNow)
        tNow = snap_t(id);
    end

    nCell{k} = nNow;
    tCell(k) = tNow;

    fprintf('  requested t = %.4f, selected t = %.4f, file = %s\n', ...
        desiredTimes(k), tNow, files(id).name);
end

gammaTag = '40';

% ============================================================
% Figure 1: top-view density snapshots
% ============================================================

fig0 = figure;
set(fig0, 'Color', 'w');
set(fig0, 'Position', topViewFigurePosition);

for k = 1:nPanels

    Z = nCell{k};
    tnow = tCell(k);

    subplot(2, 3, k);

    % Move the whole group slightly to the left
    ax = gca;
    pos = get(ax, 'Position');
    pos(1) = pos(1) - topViewLeftShift;

    if k <= 3
        pos(2) = pos(2) + topViewRowGapShift;
    else
        pos(2) = pos(2) - topViewRowGapShift;
    end

    set(ax, 'Position', pos);

    imagesc(x, y, Z);
    set(gca, 'YDir', 'normal');
    hold on;

    contour(x, y, Z, [nLevel nLevel], ...
        'k-', 'LineWidth', contourLineWidth);

    axis image;
    xlim(xwin2D);
    ylim(ywin2D);
    caxis(densityColorLimits);

    set(gca, ...
        'XTick', xtickVals2D, ...
        'YTick', ytickVals2D);

    xlabel(xLabelText, 'Interpreter', 'latex', 'FontSize', axisLabelFontSize);
    ylabel(yLabelText, 'Interpreter', 'latex', 'FontSize', axisLabelFontSize);

    title(sprintf(timeTitleFormat, tnow), ...
        'Interpreter', 'latex', 'FontSize', titleFontSize);

    set(gca, ...
        'FontSize', tickFontSize, ...
        'LineWidth', axisLineWidth, ...
        'TickLabelInterpreter', 'latex');

    box on;
end

colormap(parula);

cb = colorbar('Position', colorbarPosition);
cb.Label.String = colorbarLabelText;
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = colorbarLabelFontSize;
cb.TickLabelInterpreter = 'latex';
cb.Ticks = colorbarTicks;
cb.FontSize = colorbarFontSize;

% sgtitle('Top view of the annular density evolution', ...
%     'Interpreter', 'latex', ...
%     'FontSize', superTitleFontSize);

outfile0 = fullfile(outdir, ...
    sprintf('Figure_annulus_density_topview_gamma%s.png', gammaTag));
print(fig0, outfile0, '-dpng', sprintf('-r%d', saveResolution));

outfile0_eps = fullfile(outdir, ...
    sprintf('Figure_annulus_density_topview_gamma%s.eps', gammaTag));
print(fig0, outfile0_eps, '-depsc', sprintf('-r%d', saveResolution));

outfile0_fig = fullfile(outdir, ...
    sprintf('Figure_annulus_density_topview_gamma%s.fig', gammaTag));
savefig(fig0, outfile0_fig);

fprintf('\nSaved density top-view figure to %s\n', outfile0);

% ============================================================
% Figure 2: density cross-sections along y = 0
% ============================================================

fig1 = figure;
set(fig1, 'Color', 'w');
set(fig1, 'Position', sectionFigurePosition);

for k = 1:nPanels

    tnow = tCell(k);
    nline = nCell{k}(iy0, :);

    subplot(2, 3, k);

    % Enlarge horizontal spacing between columns
    ax = gca;
    pos = get(ax, 'Position');

    col = mod(k-1, 3) + 1;

    pos(1) = pos(1) + sectionColumnShift(col) + sectionLeftShift;
    pos(3) = sectionWidthScale * pos(3);

    set(ax, 'Position', pos);

    hold on;

    plot(x, nline, 'b-', 'LineWidth', profileLineWidth);

    xlim(xwinProfile);
    ylim(nylim);

    set(gca, ...
        'XTick', xtickValsProfile, ...
        'YTick', ytickValsProfile);

    yl = ylim;
    xcross = find_threshold_crossings(x, nline, nLevel);

    for j = 1:length(xcross)
        plot([xcross(j), xcross(j)], yl, '--', ...
            'Color', [0.4, 0.4, 0.4], ...
            'LineWidth', thresholdLineWidth);
    end

    plot(xwinProfile, [0, 0], 'k-', 'LineWidth', zeroLineWidth);

    grid on;
    box on;

    xlabel(xLabelText, 'Interpreter', 'latex', 'FontSize', axisLabelFontSize);
    ylabel(densitySectionLabelText, 'Interpreter', 'latex', ...
        'FontSize', axisLabelFontSize);

    title(sprintf(timeTitleFormat, tnow), ...
        'Interpreter', 'latex', 'FontSize', titleFontSize);

    set(gca, ...
        'FontSize', tickFontSize, ...
        'LineWidth', axisLineWidth, ...
        'TickLabelInterpreter', 'latex');
end

% sgtitle('Density cross-sections along $y=0$', ...
%     'Interpreter', 'latex', ...
%     'FontSize', superTitleFontSize);

outfile1 = fullfile(outdir, ...
    sprintf('Figure_annulus_density_sections_gamma%s.png', gammaTag));
print(fig1, outfile1, '-dpng', sprintf('-r%d', saveResolution));

outfile1_eps = fullfile(outdir, ...
    sprintf('Figure_annulus_density_sections_gamma%s.eps', gammaTag));
print(fig1, outfile1_eps, '-depsc', sprintf('-r%d', saveResolution));

outfile1_fig = fullfile(outdir, ...
    sprintf('Figure_annulus_density_sections_gamma%s.fig', gammaTag));
savefig(fig1, outfile1_fig);

fprintf('Saved density-section figure to %s\n', outfile1);

fprintf('\nAnnulus plotting finished.\n');


% ========================================================================
% Local functions
% ========================================================================

function t = parse_time_from_filename(fname)
%PARSE_TIME_FROM_FILENAME Parse time from filename.
%
% Example:
%   snap_0051_t_0p500.mat  -> 0.500
%   snap_0051_t_0p5000.mat -> 0.5000

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
%GET_SNAPSHOT_VARS Read variables from a snapshot structure.

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


function xcross = find_threshold_crossings(x, y, level)
%FIND_THRESHOLD_CROSSINGS Find approximate x-locations where y crosses level.

    xcross = [];

    for i = 1:length(x)-1

        y1 = y(i)   - level;
        y2 = y(i+1) - level;

        if y1 == 0
            xcross(end+1) = x(i); %#ok<AGROW>
        elseif y1 * y2 < 0
            xc = x(i) - y1 * (x(i+1)-x(i)) / (y2-y1);
            xcross(end+1) = xc; %#ok<AGROW>
        end
    end
end
