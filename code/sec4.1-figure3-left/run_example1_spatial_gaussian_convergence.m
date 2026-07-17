function run_example1_spatial_gaussian_convergence()
% run_example1_spatial_gaussian_convergence.m
%
% Spatial convergence test with a smooth positive Gaussian-type initial datum.
%
% Purpose:
%   Test the spatial error of the numerical scheme without using the
%   compactly supported Barenblatt profile.  A smooth strictly positive
%   Gaussian datum is used.  The reference solution is computed on the
%   finest mesh, with the same small time step as the coarser meshes.
%
% Data storage:
%   Each parameter group and each mesh are stored in separate folders.
%   Snapshots at several output times are saved as separate MAT files.
%
% Required functions on MATLAB path:
%   tumor_step_newton_n.m

clc; close all;

%% ===================== User parameters =====================

outdir = 'example1_spatial_gaussian_data';
ensure_dir(outdir);

params.a = -5.0;
params.b =  5.0;
params.T =  1.;

% Output times.  They do not have to coincide exactly with time steps.
% The nearest time step will be used.
params.outputTimes = 0:0.01:params.T;

% Smooth positive Gaussian-type initial datum:
%   n0(x) = n_bg + amp*exp(-(x-x0)^2/(2*sigma^2)).
params.n_bg  = 0.;
params.amp   = 0.4;
params.x0    = 0.0;
params.sigma = .5;

% No-growth case for the pure diffusion test.
pH   = 1.0;
Gfun = @(p) 0*p;
GpH  = 0;

% Gamma values.
params.gammaList = [3, 10, 20, 40, 100];

% Meshes for the spatial convergence test.  The finest mesh is used as
% the reference solution and is not included as a data point in the order.
params.hList = [1/4, 1/8, 1/16, 1/32, 1/64, 1/128, 1/256, 1/512, 1/1024, 1/2048];
params.hRef  = 1/2048;

% A common small time step is used for all meshes.  This makes the
% comparison mainly reflect spatial errors for the fixed time-discrete
% problem.  If the run is too slow, try 2e-4 first.
params.tau = 1.0e-5; 

% Live plot of n(x,t) during the computation.  Turn this off for batch runs.
params.showLiveNPlot = true;
params.livePlotEverySteps = 200;

% Newton options.
opts.tol = 1e-8;
opts.maxit = 100;
opts.lineSearchMaxit = 40;
opts.verbose = false;

fprintf('============================================================\n');
fprintf('Spatial convergence test with smooth Gaussian initial data\n');
fprintf('Domain = [%.3g, %.3g], T = %.6g, tau = %.4e\n', ...
    params.a, params.b, params.T, params.tau);
fprintf('Output folder: %s\n', outdir);
fprintf('============================================================\n');

%% ===================== Main loop over gamma =====================

for ig = 1:numel(params.gammaList)

    gamma = params.gammaList(ig);

    gammaDir = fullfile(outdir, sprintf('gamma_%s', num_tag(gamma)));
    ensure_dir(gammaDir);

    fprintf('\n============================================================\n');
    fprintf('Running spatial convergence test: gamma = %g\n', gamma);
    fprintf('============================================================\n');

    resultsCell = cell(numel(params.hList), 1);

    %% ---------- Run all meshes ----------
    for ih = 1:numel(params.hList)

        hTarget = params.hList(ih);

        % Use a temporary grid to determine actual N and actual h.
        Ntmp = round((params.b - params.a) / hTarget);
        hActual = (params.b - params.a) / Ntmp;

        caseTag = sprintf('gamma_%s_N_%d_h_%s_tau_%s_T_%s', ...
            num_tag(gamma), Ntmp, num_tag(hActual), ...
            num_tag(params.tau), num_tag(params.T));

        caseDir = fullfile(gammaDir, caseTag);
        ensure_dir(caseDir);
        ensure_dir(fullfile(caseDir, 'snapshots'));

        fprintf('\n------------------------------------------------------------\n');
        fprintf('Case: gamma = %g, N = %d, h = %.8g, tau = %.4e\n', ...
            gamma, Ntmp, hActual, params.tau);
        fprintf('Case folder: %s\n', caseDir);
        fprintf('------------------------------------------------------------\n');

        try
            R = run_one_gaussian_case(params, opts, gamma, pH, Gfun, GpH, ...
                hTarget, caseDir);
        catch ME
            fprintf('  This case failed: %s\n', ME.message);
            print_error_stack(ME);

            R = make_failed_result(params, gamma, hTarget, params.tau, ME.message);
        end

        resultsCell{ih} = R;

        caseFile = fullfile(caseDir, 'case_result.mat');
        save(caseFile, 'R', 'params', 'opts', 'gamma', 'pH');

        fprintf('Saved case result to:\n  %s\n', caseFile);
    end

    %% ---------- Build struct array ----------
    results = cell_to_struct_array(resultsCell);

    groupFile = fullfile(gammaDir, 'all_mesh_results.mat');
    save(groupFile, 'params', 'opts', 'gamma', 'results');

    fprintf('\nAll mesh results saved to:\n  %s\n', groupFile);

    %% ---------- Find reference solution ----------
    hValues = [results.h];
    successFlags = [results.success];

    validRef = find(successFlags & isfinite(hValues));
    if isempty(validRef)
        warning('No successful run for gamma = %g. Skip error computation.', gamma);
        continue;
    end

    [~, refLocalId] = min(abs(hValues(validRef) - params.hRef));
    refId = validRef(refLocalId);
    Rref = results(refId);

    fprintf('\nReference mesh for gamma = %g:\n', gamma);
    fprintf('  N_ref = %d, h_ref = %.8g\n', Rref.N, Rref.h);

    %% ---------- Compute errors against reference solution ----------
    summary = compute_spatial_errors(results, Rref, gamma);

    disp(summary);

    summaryFile = fullfile(gammaDir, 'spatial_error_summary.mat');
    save(summaryFile, 'summary', 'results', 'Rref', 'params', 'gamma');

    csvFile = fullfile(gammaDir, 'spatial_error_summary.csv');
    write_table_compatible(summary, csvFile);

    fprintf('Spatial error summary saved to:\n');
    fprintf('  %s\n', summaryFile);
    fprintf('  %s\n', csvFile);

    %% ---------- Plot quick check ----------
    make_quick_error_plot(summary, gammaDir, gamma);
end

fprintf('\nAll spatial Gaussian convergence tests finished.\n');

end

%% ========================================================================
%% One case
%% ========================================================================

function R = run_one_gaussian_case(params, opts, gamma, pH, Gfun, GpH, hTarget, caseDir)

    R = make_empty_result();

    a = params.a;
    b = params.b;
    T = params.T;
    tauTarget = params.tau;

    N = round((b - a) / hTarget);
    h = (b - a) / N;
    x = a + ((1:N)' - 0.5) * h;

    Nt = max(1, round(T / tauTarget));
    tau = T / Nt;
    time = (0:Nt)' * tau;

    n = gaussian_initial(x, params);
    p = n.^gamma;

    % Output steps.
    outputSteps = round(params.outputTimes(:) / tau);
    outputSteps = max(0, min(Nt, outputSteps));
    outputSteps = unique(outputSteps, 'stable');
    outputTimesActual = outputSteps * tau;
    nOut = numel(outputSteps);

    snap_n = zeros(N, nOut);
    snap_p = zeros(N, nOut);
    snap_files = cell(nOut, 1);
    snap_mass = zeros(nOut, 1);
    snap_energy = zeros(nOut, 1);

    newton_iter = zeros(Nt, 1);
    newton_res = zeros(Nt, 1);
    newton_converged = false(Nt, 1);

    min_n = zeros(Nt+1, 1);
    max_n = zeros(Nt+1, 1);
    mass = zeros(Nt+1, 1);
    energy = zeros(Nt+1, 1);

    min_n(1) = min(n);
    max_n(1) = max(n);
    mass(1) = h * sum(n);
    energy(1) = h * sum(n.^(gamma+1) / (gamma+1));

    % Save initial snapshot if requested.
    snapCounter = 0;
    if any(outputSteps == 0)
        snapCounter = snapCounter + 1;
        [snap_n, snap_p, snap_files, snap_mass, snap_energy] = save_snapshot( ...
            snapCounter, x, n, p, gamma, h, tau, T, 0, 0, ...
            outputTimesActual(snapCounter), caseDir, snap_n, snap_p, ...
            snap_files, snap_mass, snap_energy);
    end

    tic;
    reportEvery = max(1, round(Nt / 10));
    showLiveNPlot = get_param(params, 'showLiveNPlot', false);
    livePlotEverySteps = max(1, round(get_param(params, ...
        'livePlotEverySteps', reportEvery)));

    if showLiveNPlot
        livePlot = initialize_live_n_plot(x, n, gamma, N, h, 0);
    else
        livePlot = [];
    end

    for k = 1:Nt

        [n, p, ~, info] = tumor_step_newton_n(n, h, tau, gamma, pH, Gfun, GpH, opts);
        n = n(:);
        p = p(:);

        if isfield(info, 'converged')
            newton_converged(k) = info.converged;
        else
            newton_converged(k) = true;
        end

        if isfield(info, 'iter')
            newton_iter(k) = info.iter;
        else
            newton_iter(k) = NaN;
        end

        if isfield(info, 'residual')
            newton_res(k) = info.residual;
        else
            newton_res(k) = NaN;
        end

        if isfield(info, 'converged') && ~info.converged
            fprintf('    Newton not fully converged at step %d/%d, residual = %.4e\n', ...
                k, Nt, newton_res(k));
        end

        min_n(k+1) = min(n);
        max_n(k+1) = max(n);
        mass(k+1) = h * sum(n);
        energy(k+1) = h * sum(n.^(gamma+1) / (gamma+1));

        if any(outputSteps == k)
            snapCounter = snapCounter + 1;
            [snap_n, snap_p, snap_files, snap_mass, snap_energy] = save_snapshot( ...
                snapCounter, x, n, p, gamma, h, tau, T, k, k*tau, ...
                outputTimesActual(snapCounter), caseDir, snap_n, snap_p, ...
                snap_files, snap_mass, snap_energy);
        end

        if showLiveNPlot && (k == 1 || mod(k, livePlotEverySteps) == 0 || k == Nt)
            update_live_n_plot(livePlot, n, k*tau, k, Nt, ...
                newton_iter(k), newton_res(k));
        end

        if mod(k, reportEvery) == 0 || k == Nt
            fprintf('  step %6d/%6d, t = %.5f, Newton it = %d, res = %.3e\n', ...
                k, Nt, k*tau, newton_iter(k), newton_res(k));
        end
    end

    runtime = toc;

    R.success = true;
    R.error_message = '';
    R.a = a;
    R.b = b;
    R.N = N;
    R.h = h;
    R.x = x;
    R.tau = tau;
    R.tau_target = tauTarget;
    R.T = T;
    R.Nt = Nt;
    R.gamma = gamma;
    R.initial_type = 'smooth_positive_gaussian';
    R.initial_params = struct('n_bg', params.n_bg, 'amp', params.amp, ...
        'x0', params.x0, 'sigma', params.sigma);
    R.output_steps = outputSteps;
    R.output_times = outputTimesActual;
    R.snap_n = snap_n;
    R.snap_p = snap_p;
    R.snap_files = snap_files;
    R.snap_mass = snap_mass;
    R.snap_energy = snap_energy;
    R.final_n = n;
    R.final_p = p;
    R.time = time;
    R.mass = mass;
    R.energy = energy;
    R.min_n = min_n;
    R.max_n = max_n;
    R.newton_iter = newton_iter;
    R.newton_res = newton_res;
    R.newton_converged = newton_converged;
    R.runtime = runtime;
end

%% ========================================================================
%% Snapshot helper
%% ========================================================================

function [snap_n, snap_p, snap_files, snap_mass, snap_energy] = save_snapshot( ...
    idx, x, n, p, gamma, h, tau, T, step, t, outputTime, caseDir, ...
    snap_n, snap_p, snap_files, snap_mass, snap_energy)

    snap_n(:, idx) = n;
    snap_p(:, idx) = p;
    snap_mass(idx) = h * sum(n);
    snap_energy(idx) = h * sum(n.^(gamma+1) / (gamma+1));

    N = numel(x);
    snapDir = fullfile(caseDir, 'snapshots');
    ensure_dir(snapDir);

    snapFile = fullfile(snapDir, sprintf( ...
        'snapshot_gamma_%s_N_%d_h_%s_tau_%s_t_%s.mat', ...
        num_tag(gamma), N, num_tag(h), num_tag(tau), num_tag(t)));

    save(snapFile, 'x', 'n', 'p', 'gamma', 'N', 'h', 'tau', 'T', ...
        'step', 't', 'outputTime');

    snap_files{idx} = snapFile;
end

%% ========================================================================
%% Error computation
%% ========================================================================

function summary = compute_spatial_errors(results, Rref, gamma)

    nCases = numel(results);

    h_col = [];
    N_col = [];
    tau_col = [];
    err_L1_T_col = [];
    err_Linf_T_col = [];
    err_L1_max_col = [];
    err_Linf_max_col = [];
    order_L1_T_col = [];
    order_Linf_T_col = [];
    order_L1_max_col = [];
    order_Linf_max_col = [];

    for i = 1:nCases

        R = results(i);

        if ~R.success
            continue;
        end

        % The reference mesh itself is not included in the error table.
        if abs(R.h - Rref.h) <= 1e-12 * max(1, abs(Rref.h))
            continue;
        end

        if numel(R.output_times) ~= numel(Rref.output_times)
            error('Output time arrays have different lengths.');
        end

        nOut = numel(R.output_times);
        errL1_time = zeros(nOut, 1);
        errLinf_time = zeros(nOut, 1);

        for j = 1:nOut
            nRefInterp = interp1(Rref.x, Rref.snap_n(:, j), R.x, 'pchip');
            nRefInterp = nRefInterp(:);

            diffj = R.snap_n(:, j) - nRefInterp;
            errL1_time(j) = R.h * sum(abs(diffj));
            errLinf_time(j) = max(abs(diffj));
        end

        h_col(end+1,1) = R.h; %#ok<AGROW>
        N_col(end+1,1) = R.N; %#ok<AGROW>
        tau_col(end+1,1) = R.tau; %#ok<AGROW>
        err_L1_T_col(end+1,1) = errL1_time(end); %#ok<AGROW>
        err_Linf_T_col(end+1,1) = errLinf_time(end); %#ok<AGROW>
        err_L1_max_col(end+1,1) = max(errL1_time); %#ok<AGROW>
        err_Linf_max_col(end+1,1) = max(errLinf_time); %#ok<AGROW>
    end

    % Sort from coarse to fine.
    [h_col, idx] = sort(h_col, 'descend');
    N_col = N_col(idx);
    tau_col = tau_col(idx);
    err_L1_T_col = err_L1_T_col(idx);
    err_Linf_T_col = err_Linf_T_col(idx);
    err_L1_max_col = err_L1_max_col(idx);
    err_Linf_max_col = err_Linf_max_col(idx);

    order_L1_T_col = compute_orders_h(h_col, err_L1_T_col);
    order_Linf_T_col = compute_orders_h(h_col, err_Linf_T_col);
    order_L1_max_col = compute_orders_h(h_col, err_L1_max_col);
    order_Linf_max_col = compute_orders_h(h_col, err_Linf_max_col);

    gamma_col = gamma * ones(numel(h_col), 1);
    h_ref_col = Rref.h * ones(numel(h_col), 1);
    N_ref_col = Rref.N * ones(numel(h_col), 1);

    summary = table(gamma_col, N_col, h_col, tau_col, N_ref_col, h_ref_col, ...
        err_L1_T_col, order_L1_T_col, err_Linf_T_col, order_Linf_T_col, ...
        err_L1_max_col, order_L1_max_col, err_Linf_max_col, order_Linf_max_col, ...
        'VariableNames', {'gamma', 'N', 'h', 'tau', 'N_ref', 'h_ref', ...
        'err_L1_T', 'order_L1_T', 'err_Linf_T', 'order_Linf_T', ...
        'err_L1_max', 'order_L1_max', 'err_Linf_max', 'order_Linf_max'});
end

function orders = compute_orders_h(h, err)

    orders = NaN(size(err));

    for j = 2:numel(err)
        if err(j) > 0 && err(j-1) > 0 && h(j) > 0 && h(j-1) > 0
            orders(j) = log(err(j-1) / err(j)) / log(h(j-1) / h(j));
        end
    end
end

%% ========================================================================
%% Initial condition and result structures
%% ========================================================================

function n0 = gaussian_initial(x, params)

    n0 = params.n_bg + params.amp * exp(-((x - params.x0).^2) / (2 * params.sigma^2));
    n0 = n0(:);
end

function livePlot = initialize_live_n_plot(x, n, gamma, N, h, t)

    livePlot.fig = figure('Name', sprintf('Live n, gamma=%g, N=%d', gamma, N));
    set(livePlot.fig, 'Color', 'w');

    livePlot.ax = axes('Parent', livePlot.fig);
    livePlot.line = plot(livePlot.ax, x, n, 'b-', 'LineWidth', 2.0);

    xlabel(livePlot.ax, '$x$', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel(livePlot.ax, '$n(x,t)$', 'Interpreter', 'latex', 'FontSize', 16);
    title(livePlot.ax, sprintf('$\\gamma=%g,\\ N=%d,\\ h=%.3g,\\ t=%.5f$', ...
        gamma, N, h, t), 'Interpreter', 'latex', 'FontSize', 16);

    set(livePlot.ax, 'FontSize', 14, 'LineWidth', 1.0, ...
        'TickLabelInterpreter', 'latex');
    grid(livePlot.ax, 'on');
    box(livePlot.ax, 'on');
    xlim(livePlot.ax, [x(1), x(end)]);
    ylim(livePlot.ax, [0, max(1.0e-12, 1.05 * max(n))]);

    livePlot.gamma = gamma;
    livePlot.N = N;
    livePlot.h = h;

    drawnow;
end

function update_live_n_plot(livePlot, n, t, step, Nt, newtonIter, newtonRes)

    if isempty(livePlot) || ~ishandle(livePlot.fig) || ~ishandle(livePlot.line)
        return;
    end

    set(livePlot.line, 'YData', n);

    ymax = max(n);
    if ymax <= 0
        ymax = 1.0e-12;
    end
    ylim(livePlot.ax, [0, 1.05 * ymax]);

    title(livePlot.ax, sprintf([ ...
        '$\\gamma=%g,\\ N=%d,\\ h=%.3g,\\ t=%.5f,\\ ', ...
        'step=%d/%d,\\ Newton=%g,\\ res=%.2e$'], ...
        livePlot.gamma, livePlot.N, livePlot.h, t, ...
        step, Nt, newtonIter, newtonRes), ...
        'Interpreter', 'latex', 'FontSize', 16);

    drawnow limitrate;
end

function R = make_empty_result()

    R.success = false;
    R.error_message = '';
    R.a = NaN;
    R.b = NaN;
    R.N = NaN;
    R.h = NaN;
    R.x = [];
    R.tau = NaN;
    R.tau_target = NaN;
    R.T = NaN;
    R.Nt = NaN;
    R.gamma = NaN;
    R.initial_type = '';
    R.initial_params = [];
    R.output_steps = [];
    R.output_times = [];
    R.snap_n = [];
    R.snap_p = [];
    R.snap_files = [];
    R.snap_mass = [];
    R.snap_energy = [];
    R.final_n = [];
    R.final_p = [];
    R.time = [];
    R.mass = [];
    R.energy = [];
    R.min_n = [];
    R.max_n = [];
    R.newton_iter = [];
    R.newton_res = [];
    R.newton_converged = [];
    R.runtime = NaN;
end

function R = make_failed_result(params, gamma, h, tau, msg)

    R = make_empty_result();
    R.success = false;
    R.error_message = msg;
    R.a = params.a;
    R.b = params.b;
    R.h = h;
    R.tau = tau;
    R.T = params.T;
    R.gamma = gamma;
end

function S = cell_to_struct_array(C)

    C = C(:);
    keep = true(size(C));

    for i = 1:numel(C)
        if isempty(C{i})
            keep(i) = false;
        end
    end

    C = C(keep);

    if isempty(C)
        S = make_empty_result();
        S = S([]);
        return;
    end

    template = make_empty_result();
    allFields = fieldnames(template);

    for i = 1:numel(C)
        fi = fieldnames(C{i});
        for j = 1:numel(fi)
            if ~ismember(fi{j}, allFields)
                allFields{end+1,1} = fi{j}; %#ok<AGROW>
                template.(fi{j}) = [];
            end
        end
    end

    S = repmat(template, 1, numel(C));

    for i = 1:numel(C)
        Si = template;
        fi = fieldnames(C{i});
        for j = 1:numel(fi)
            Si.(fi{j}) = C{i}.(fi{j});
        end
        S(i) = Si;
    end
end

%% ========================================================================
%% Plot and utilities
%% ========================================================================

function make_quick_error_plot(summary, gammaDir, gamma)

    if isempty(summary)
        return;
    end

    fig = figure('Color', 'w');
    hold on; box on;

    loglog(summary.h, summary.err_L1_T, 'o-', ...
        'LineWidth', 1.5, 'MarkerSize', 7, ...
        'DisplayName', '$L^1$ error at $T$');

    loglog(summary.h, summary.err_Linf_T, 's-', ...
        'LineWidth', 1.5, 'MarkerSize', 7, ...
        'DisplayName', '$L^\infty$ error at $T$');

    href = summary.h;
    ref = summary.err_L1_T(end) * (href / href(end)).^2;
    loglog(href, ref, 'k--', 'LineWidth', 1.2, 'DisplayName', '$O(h^2)$ reference');

    xlabel('$h$', 'Interpreter', 'latex');
    ylabel('error against finest-grid reference', 'Interpreter', 'latex');
    title(sprintf('Spatial convergence, $\\gamma=%g$', gamma), 'Interpreter', 'latex');
    legend('Location', 'best', 'Interpreter', 'latex');
    grid on;
    set(gca, 'FontSize', 12);

    figFile = fullfile(gammaDir, sprintf('spatial_error_gamma_%s.png', num_tag(gamma)));
    print(fig, figFile, '-dpng', '-r300');
end

function ensure_dir(folder)

    if ~exist(folder, 'dir')
        mkdir(folder);
    end
end

function value = get_param(params, fieldName, defaultValue)

    if isfield(params, fieldName)
        value = params.(fieldName);
    else
        value = defaultValue;
    end
end

function write_table_compatible(T, filename)

    try
        writetable(T, filename);
    catch
        fprintf('writetable failed. Saving table as MAT file instead.\n');
        [folder, name, ~] = fileparts(filename);
        save(fullfile(folder, [name, '.mat']), 'T');
    end
end

function print_error_stack(ME)

    try
        st = ME.stack;
    catch
        st = [];
    end

    if ~isempty(st)
        fprintf('  Error stack:\n');
        for k = 1:numel(st)
            fprintf('    %s, line %d\n', st(k).name, st(k).line);
        end
    end
end

function s = num_tag(x)

    s = sprintf('%.12g', x);
    s = strrep(s, '.', 'p');
    s = strrep(s, '-', 'm');
    s = strrep(s, '+', '');
end
