function run_barenblatt_configurable()
%RUN_BARENBLATT_CONFIGURABLE Configurable 1D Barenblatt / growth tests.
%
% This driver replaces the old Gaussian-initial-data convergence driver.
% It keeps the nonlinear step in tumor_step_newton_n.m unchanged: the solver
% still uses n as the Newton unknown, p=n^gamma, and a damped Newton line
% search that preserves nonnegativity.
%
% Required function on MATLAB path:
%   tumor_step_newton_n.m
%
% Main choices implemented here:
%   1. Initial datum is the time-shifted Barenblatt profile n_ex(x,0).
%   2. gamma can be 3, 20, or both.
%   3. C is chosen automatically by default: C=1 for gamma=3, C=0.1 for gamma=20.
%   4. tauList may contain one value or several values.
%   5. NList may contain one value or several values.
%   6. GCaseList may be {'zero'}, {'linear'}, or {'zero','linear'}, where
%          'zero'   means G(p)=0,
%          'linear' means G(p)=1-p, pH=1.
%
% Error convention:
%   * For G(p)=0, the Barenblatt solution is exact, so errors are computed
%     against n_ex(x,T).
%   * For G(p)=1-p, the Barenblatt profile is used only as the initial data.
%     There is no Barenblatt exact solution for this growth equation, so the
%     temporal/spatial summaries use the finest available numerical run as
%     reference inside each fixed-N or fixed-tau group.

clc; close all;

%% ===================== User parameters =====================

outdir = 'barenblatt_configurable_data';
ensure_dir(outdir);

params = struct();

% Domain and final time.  The paper uses [-5,5] for the 1D tests.
params.a = -5.0;
params.b =  5.0;
params.T =  1;

% Barenblatt time shift.
params.t0 = 2;

% Output times.  They do not have to coincide exactly with time steps;
% the nearest time step is used.  The code always also includes t=0 and T.
params.outputTimes = linspace(0, params.T, 101);

% gamma choices: use [3], [20], or [3 20].
params.gammaList = [3 20];

% Barenblatt constants.  These defaults reproduce the paper benchmark.
params.CForGamma3  = 1.0;
params.CForGamma20 = 0.1;

% Growth choices:
%   {'zero'}          -> only G(p)=0
%   {'linear'}        -> only G(p)=1-p
%   {'zero','linear'} -> run both
params.GCaseList = {'zero', 'linear'};

% Spatial grids are specified by the number of cell centers N.
% On [-5,5], N=640 gives h=1/64, as in the paper profile test.
% For spatial refinement, for example use: [160 320 640 1280].
params.NList = 4096;

% Time steps.  For the paper profile test with N=640, h=1/64 and tau=0.01*h.
% For temporal refinement, for example use: [1e-3 5e-4 2.5e-4 1.25e-4].
h_paper = (params.b - params.a) / 640;
params.tauList = 0.1 * [1, 1/2, 1/4, 1/8, 1/16, 1/32, 1/64, 1/128, 1/256, 1/512];

% Diagnostics and plotting.
params.saveSnapshots = true;
params.makeFinalProfilePlot = true;
params.showLiveNPlot = false;
params.livePlotEverySteps = 100;

% Newton options.  Keep these consistent with Appendix A.
opts = struct();
opts.tol = 1e-10;
opts.maxit = 100;
opts.lineSearchMaxit = 40;
opts.tolH = 1e-12;
opts.armijo = 1e-4;
opts.verbose = false;

%% ===================== Normalize user lists =====================

params.gammaList = params.gammaList(:)';
params.NList = params.NList(:)';
params.tauList = params.tauList(:)';
params.GCaseList = normalize_cellstr(params.GCaseList);
params.outputTimes = sanitize_output_times(params.outputTimes, params.T);

validate_user_parameters(params);

fprintf('============================================================\n');
fprintf('Configurable Barenblatt benchmark / growth driver\n');
fprintf('Domain = [%.3g, %.3g], T = %.6g, t0 = %.6g\n', ...
    params.a, params.b, params.T, params.t0);
fprintf('gammaList = %s\n', mat2str(params.gammaList));
fprintf('NList     = %s\n', mat2str(params.NList));
fprintf('tauList   = %s\n', mat2str(params.tauList));
fprintf('GCaseList = %s\n', strjoin_compat(params.GCaseList, ', '));
fprintf('Output folder: %s\n', outdir);
fprintf('============================================================\n');

%% ===================== Main loops =====================

for iG = 1:numel(params.GCaseList)

    GCase = params.GCaseList{iG};
    [Gfun, pH, GpH, GTag, GLabel, exactAvailable, energyPH] = get_growth_case(GCase);

    for ig = 1:numel(params.gammaList)

        gamma = params.gammaList(ig);
        C = get_barenblatt_C(params, gamma);

        groupDir = fullfile(outdir, GTag, ...
            sprintf('gamma_%s_C_%s', num_tag(gamma), num_tag(C)));
        ensure_dir(groupDir);

        fprintf('\n============================================================\n');
        fprintf('Running group: %s, gamma = %g, C = %g\n', GLabel, gamma, C);
        if exactAvailable
            fprintf('Error mode: exact Barenblatt profile.\n');
        else
            fprintf('Error mode: numerical reference; Barenblatt is only the initial datum.\n');
        end
        fprintf('Group folder: %s\n', groupDir);
        fprintf('============================================================\n');

        nCases = numel(params.NList) * numel(params.tauList);
        resultsCell = cell(nCases, 1);
        icase = 0;

        for iN = 1:numel(params.NList)
            NTarget = params.NList(iN);

            for itau = 1:numel(params.tauList)
                tauTarget = params.tauList(itau);

                hActual = (params.b - params.a) / NTarget;
                NtActual = max(1, round(params.T / tauTarget));
                tauActual = params.T / NtActual;

                caseTag = sprintf('%s_gamma_%s_C_%s_N_%d_h_%s_tau_%s_T_%s', ...
                    GTag, num_tag(gamma), num_tag(C), NTarget, ...
                    num_tag(hActual), num_tag(tauActual), num_tag(params.T));
                caseDir = fullfile(groupDir, caseTag);
                ensure_dir(caseDir);
                if params.saveSnapshots
                    ensure_dir(fullfile(caseDir, 'snapshots'));
                end

                fprintf('\n------------------------------------------------------------\n');
                fprintf('Case %d/%d: %s, gamma = %g, C = %g, N = %d, h = %.8g, tau = %.4e, Nt = %d\n', ...
                    icase + 1, nCases, GLabel, gamma, C, NTarget, hActual, tauActual, NtActual);
                fprintf('Case folder: %s\n', caseDir);
                fprintf('------------------------------------------------------------\n');

                icase = icase + 1;

                try
                    R = run_one_barenblatt_case(params, opts, gamma, C, ...
                        GCase, GTag, GLabel, pH, Gfun, GpH, exactAvailable, ...
                        energyPH, NTarget, tauTarget, caseDir);
                catch ME
                    fprintf('  This case failed: %s\n', ME.message);
                    print_error_stack(ME);
                    R = make_failed_result(params, gamma, C, GCase, GTag, GLabel, ...
                        pH, energyPH, exactAvailable, NTarget, tauTarget, ME.message);
                end

                resultsCell{icase} = R;

                caseFile = fullfile(caseDir, 'case_result.mat');
                save(caseFile, 'R', 'params', 'opts', 'gamma', 'C', ...
                    'GCase', 'GTag', 'GLabel', 'pH', 'GpH', 'energyPH');
                fprintf('Saved case result to:\n  %s\n', caseFile);
            end
        end

        results = cell_to_struct_array(resultsCell);
        allFile = fullfile(groupDir, 'all_case_results.mat');
        save(allFile, 'params', 'opts', 'gamma', 'C', 'GCase', 'GTag', ...
            'GLabel', 'pH', 'GpH', 'energyPH', 'exactAvailable', 'results');
        fprintf('\nAll case results saved to:\n  %s\n', allFile);

        caseTable = build_case_table(results);
        caseCsv = fullfile(groupDir, 'case_table.csv');
        write_table_compatible(caseTable, caseCsv);
        save(fullfile(groupDir, 'case_table.mat'), 'caseTable');

        temporalSummary = build_temporal_summary(results, exactAvailable);
        spatialSummary = build_spatial_summary(results, exactAvailable);

        temporalCsv = fullfile(groupDir, 'temporal_error_summary.csv');
        spatialCsv  = fullfile(groupDir, 'spatial_error_summary.csv');
        write_table_compatible(temporalSummary, temporalCsv);
        write_table_compatible(spatialSummary, spatialCsv);
        save(fullfile(groupDir, 'error_summaries.mat'), ...
            'temporalSummary', 'spatialSummary', 'caseTable', 'results');

        fprintf('\nCase table and error summaries saved to:\n');
        fprintf('  %s\n', caseCsv);
        fprintf('  %s\n', temporalCsv);
        fprintf('  %s\n', spatialCsv);

        disp('Temporal summary:');
        disp(temporalSummary);
        disp('Spatial summary:');
        disp(spatialSummary);

        if params.makeFinalProfilePlot
            make_group_profile_plot(results, groupDir, exactAvailable);
            make_error_summary_plots(temporalSummary, spatialSummary, groupDir, exactAvailable);
        end
    end
end

fprintf('\nAll configurable Barenblatt runs finished.\n');

end

%% ========================================================================
%% One case
%% ========================================================================

function R = run_one_barenblatt_case(params, opts, gamma, C, ...
    GCase, GTag, GLabel, pH, Gfun, GpH, exactAvailable, energyPH, ...
    NTarget, tauTarget, caseDir)

    R = make_empty_result();

    a = params.a;
    b = params.b;
    T = params.T;
    t0 = params.t0;

    N = round(NTarget);
    h = (b - a) / N;
    x = a + ((1:N)' - 0.5) * h;

    Nt = max(1, round(T / tauTarget));
    tau = T / Nt;
    time = (0:Nt)' * tau;

    n = barenblatt_exact(x, 0.0, gamma, C, t0);
    n = max(n(:), 0.0);
    p = n.^gamma;

    if strcmpi(GTag, 'G1mp') && max(p) > pH * (1 + 1e-12)
        warning(['For G(p)=1-p, the homeostatic upper-bound assumption ', ...
            'p0 <= pH is not satisfied by this Barenblatt initial datum ', ...
            '(max p0 = %.6g, pH = %.6g). The run is still performed, ', ...
            'but the Barenblatt profile is only used as initial data.'], max(p), pH);
    end

    outputSteps = round(params.outputTimes(:) / tau);
    outputSteps = max(0, min(Nt, outputSteps));
    outputSteps = unique(outputSteps, 'stable');
    outputTimesActual = outputSteps * tau;
    nOut = numel(outputSteps);

    snap_n = zeros(N, nOut);
    snap_p = zeros(N, nOut);
    snap_n_exact = NaN(N, nOut);
    snap_files = cell(nOut, 1);
    snap_mass = zeros(nOut, 1);
    snap_energy = zeros(nOut, 1);
    snap_L1_exact = NaN(nOut, 1);
    snap_Linf_exact = NaN(nOut, 1);

    newton_iter = zeros(Nt, 1);
    newton_res = zeros(Nt, 1);
    newton_converged = false(Nt, 1);

    min_n = zeros(Nt+1, 1);
    max_n = zeros(Nt+1, 1);
    min_p = zeros(Nt+1, 1);
    max_p = zeros(Nt+1, 1);
    mass = zeros(Nt+1, 1);
    energy = zeros(Nt+1, 1);

    min_n(1) = min(n);
    max_n(1) = max(n);
    min_p(1) = min(p);
    max_p(1) = max(p);
    mass(1) = h * sum(n);
    energy(1) = discrete_energy(n, h, gamma, energyPH);

    snapCounter = 0;
    if any(outputSteps == 0)
        snapCounter = snapCounter + 1;
        [snap_n, snap_p, snap_n_exact, snap_files, snap_mass, snap_energy, ...
            snap_L1_exact, snap_Linf_exact] = save_snapshot( ...
            snapCounter, x, n, p, gamma, C, t0, h, tau, T, 0, 0, ...
            outputTimesActual(snapCounter), caseDir, GTag, GLabel, ...
            exactAvailable, energyPH, params.saveSnapshots, ...
            snap_n, snap_p, snap_n_exact, snap_files, snap_mass, snap_energy, ...
            snap_L1_exact, snap_Linf_exact);
    end

    tic;
    reportEvery = max(1, round(Nt / 10));
    showLiveNPlot = get_param(params, 'showLiveNPlot', false);
    livePlotEverySteps = max(1, round(get_param(params, ...
        'livePlotEverySteps', reportEvery)));

    if showLiveNPlot
        livePlot = initialize_live_n_plot(x, n, gamma, C, GLabel, N, h, 0);
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
        min_p(k+1) = min(p);
        max_p(k+1) = max(p);
        mass(k+1) = h * sum(n);
        energy(k+1) = discrete_energy(n, h, gamma, energyPH);

        if any(outputSteps == k)
            snapCounter = snapCounter + 1;
            [snap_n, snap_p, snap_n_exact, snap_files, snap_mass, snap_energy, ...
                snap_L1_exact, snap_Linf_exact] = save_snapshot( ...
                snapCounter, x, n, p, gamma, C, t0, h, tau, T, k, k*tau, ...
                outputTimesActual(snapCounter), caseDir, GTag, GLabel, ...
                exactAvailable, energyPH, params.saveSnapshots, ...
                snap_n, snap_p, snap_n_exact, snap_files, snap_mass, snap_energy, ...
                snap_L1_exact, snap_Linf_exact);
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

    if exactAvailable
        nExactFinal = barenblatt_exact(x, T, gamma, C, t0);
        diffFinal = n - nExactFinal;
        errL1ExactT = h * sum(abs(diffFinal));
        errLinfExactT = max(abs(diffFinal));
        errL1ExactMax = max(snap_L1_exact);
        errLinfExactMax = max(snap_Linf_exact);
    else
        nExactFinal = NaN(size(n));
        errL1ExactT = NaN;
        errLinfExactT = NaN;
        errL1ExactMax = NaN;
        errLinfExactMax = NaN;
    end

    if abs(mass(1)) > 0
        relMassChange = (mass(end) - mass(1)) / mass(1);
    else
        relMassChange = NaN;
    end

    if abs(energy(1)) > 0
        energyRatio = energy(end) / energy(1);
    else
        energyRatio = NaN;
    end

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
    R.C = C;
    R.t0 = t0;
    R.G_case = GCase;
    R.G_tag = GTag;
    R.G_label = GLabel;
    R.pH = pH;
    R.GpH = GpH;
    R.energy_pH = energyPH;
    R.exact_available = exactAvailable;
    R.initial_type = 'barenblatt_shifted';
    R.output_steps = outputSteps;
    R.output_times = outputTimesActual;
    R.snap_n = snap_n;
    R.snap_p = snap_p;
    R.snap_n_exact = snap_n_exact;
    R.snap_files = snap_files;
    R.snap_mass = snap_mass;
    R.snap_energy = snap_energy;
    R.snap_L1_exact = snap_L1_exact;
    R.snap_Linf_exact = snap_Linf_exact;
    R.final_n = n;
    R.final_p = p;
    R.final_n_exact = nExactFinal;
    R.time = time;
    R.mass = mass;
    R.energy = energy;
    R.min_n = min_n;
    R.max_n = max_n;
    R.min_p = min_p;
    R.max_p = max_p;
    R.rel_mass_change = relMassChange;
    R.energy_ratio = energyRatio;
    R.err_L1_exact_T = errL1ExactT;
    R.err_Linf_exact_T = errLinfExactT;
    R.err_L1_exact_max = errL1ExactMax;
    R.err_Linf_exact_max = errLinfExactMax;
    R.newton_iter = newton_iter;
    R.newton_res = newton_res;
    R.newton_converged = newton_converged;
    R.runtime = runtime;

    if get_param(params, 'makeFinalProfilePlot', false)
        make_case_profile_plot(R, caseDir);
    end
end

%% ========================================================================
%% Snapshot helper
%% ========================================================================

function [snap_n, snap_p, snap_n_exact, snap_files, snap_mass, snap_energy, ...
    snap_L1_exact, snap_Linf_exact] = save_snapshot( ...
    idx, x, n, p, gamma, C, t0, h, tau, T, step, t, outputTime, caseDir, ...
    GTag, GLabel, exactAvailable, energyPH, doSave, ...
    snap_n, snap_p, snap_n_exact, snap_files, snap_mass, snap_energy, ...
    snap_L1_exact, snap_Linf_exact)

    snap_n(:, idx) = n;
    snap_p(:, idx) = p;
    snap_mass(idx) = h * sum(n);
    snap_energy(idx) = discrete_energy(n, h, gamma, energyPH);

    if exactAvailable
        nExact = barenblatt_exact(x, t, gamma, C, t0);
        snap_n_exact(:, idx) = nExact;
        diffNow = n - nExact;
        snap_L1_exact(idx) = h * sum(abs(diffNow));
        snap_Linf_exact(idx) = max(abs(diffNow));
    else
        nExact = NaN(size(n));
    end

    if ~doSave
        snap_files{idx} = '';
        return;
    end

    N = numel(x);
    snapDir = fullfile(caseDir, 'snapshots');
    ensure_dir(snapDir);

    snapFile = fullfile(snapDir, sprintf( ...
        'snapshot_%s_gamma_%s_C_%s_N_%d_h_%s_tau_%s_t_%s.mat', ...
        GTag, num_tag(gamma), num_tag(C), N, num_tag(h), num_tag(tau), num_tag(t)));

    save(snapFile, 'x', 'n', 'p', 'nExact', 'gamma', 'C', 't0', ...
        'GTag', 'GLabel', 'N', 'h', 'tau', 'T', 'step', 't', 'outputTime', ...
        'energyPH');

    snap_files{idx} = snapFile;
end

%% ========================================================================
%% Error summaries
%% ========================================================================

function temporalSummary = build_temporal_summary(results, exactAvailable)

    valid = find_success_indices(results);
    if isempty(valid)
        temporalSummary = table();
        return;
    end

    NValues = unique([results(valid).N]);

    gamma_col = [];
    C_col = [];
    N_col = [];
    h_col = [];
    tau_col = [];
    Nt_col = [];
    tau_ref_col = [];
    err_L1_T_col = [];
    order_L1_T_col = [];
    err_Linf_T_col = [];
    order_Linf_T_col = [];
    error_mode_col = {};
    G_label_col = {};

    for iN = 1:numel(NValues)
        N0 = NValues(iN);
        ids = valid([results(valid).N] == N0);
        if numel(ids) < 2
            continue;
        end

        [~, perm] = sort([results(ids).tau], 'descend');
        ids = ids(perm);

        if exactAvailable
            stepVals = [results(ids).tau]';
            errL1 = [results(ids).err_L1_exact_T]';
            errLinf = [results(ids).err_Linf_exact_T]';
            tauRef = NaN;
            modeText = 'exact Barenblatt';
        else
            [~, refLocal] = min([results(ids).tau]);
            refId = ids(refLocal);
            Rref = results(refId);
            ids(ids == refId) = [];
            if isempty(ids)
                continue;
            end
            [~, perm2] = sort([results(ids).tau], 'descend');
            ids = ids(perm2);
            stepVals = [results(ids).tau]';
            errL1 = zeros(numel(ids), 1);
            errLinf = zeros(numel(ids), 1);
            for j = 1:numel(ids)
                R = results(ids(j));
                diffj = R.final_n - Rref.final_n;
                errL1(j) = R.h * sum(abs(diffj));
                errLinf(j) = max(abs(diffj));
            end
            tauRef = Rref.tau;
            modeText = sprintf('reference tau=%g', tauRef);
        end

        orderL1 = compute_orders(stepVals, errL1);
        orderLinf = compute_orders(stepVals, errLinf);

        for j = 1:numel(ids)
            R = results(ids(j));
            gamma_col(end+1,1) = R.gamma; %#ok<AGROW>
            C_col(end+1,1) = R.C; %#ok<AGROW>
            N_col(end+1,1) = R.N; %#ok<AGROW>
            h_col(end+1,1) = R.h; %#ok<AGROW>
            tau_col(end+1,1) = R.tau; %#ok<AGROW>
            Nt_col(end+1,1) = R.Nt; %#ok<AGROW>
            tau_ref_col(end+1,1) = tauRef; %#ok<AGROW>
            err_L1_T_col(end+1,1) = errL1(j); %#ok<AGROW>
            order_L1_T_col(end+1,1) = orderL1(j); %#ok<AGROW>
            err_Linf_T_col(end+1,1) = errLinf(j); %#ok<AGROW>
            order_Linf_T_col(end+1,1) = orderLinf(j); %#ok<AGROW>
            error_mode_col{end+1,1} = modeText; %#ok<AGROW>
            G_label_col{end+1,1} = R.G_label; %#ok<AGROW>
        end
    end

    if isempty(gamma_col)
        temporalSummary = table();
        return;
    end

    temporalSummary = table(G_label_col, gamma_col, C_col, N_col, h_col, ...
        tau_col, Nt_col, tau_ref_col, err_L1_T_col, order_L1_T_col, ...
        err_Linf_T_col, order_Linf_T_col, error_mode_col, ...
        'VariableNames', {'G','gamma','C','N','h','tau','Nt','tau_ref', ...
        'err_L1_T','order_L1_T','err_Linf_T','order_Linf_T','error_mode'});
end

function spatialSummary = build_spatial_summary(results, exactAvailable)

    valid = find_success_indices(results);
    if isempty(valid)
        spatialSummary = table();
        return;
    end

    tauTargets = unique(round_sig([results(valid).tau_target], 14));

    gamma_col = [];
    C_col = [];
    N_col = [];
    h_col = [];
    tau_col = [];
    N_ref_col = [];
    h_ref_col = [];
    err_L1_T_col = [];
    order_L1_T_col = [];
    err_Linf_T_col = [];
    order_Linf_T_col = [];
    error_mode_col = {};
    G_label_col = {};

    for itau = 1:numel(tauTargets)
        tau0 = tauTargets(itau);
        ids = [];
        for k = 1:numel(valid)
            idx = valid(k);
            if round_sig(results(idx).tau_target, 14) == tau0
                ids(end+1) = idx; %#ok<AGROW>
            end
        end

        if numel(ids) < 2
            continue;
        end

        [~, perm] = sort([results(ids).h], 'descend');
        ids = ids(perm);

        if exactAvailable
            stepVals = [results(ids).h]';
            errL1 = [results(ids).err_L1_exact_T]';
            errLinf = [results(ids).err_Linf_exact_T]';
            Nref = NaN;
            href = NaN;
            modeText = 'exact Barenblatt';
        else
            [~, refLocal] = min([results(ids).h]);
            refId = ids(refLocal);
            Rref = results(refId);
            ids(ids == refId) = [];
            if isempty(ids)
                continue;
            end
            [~, perm2] = sort([results(ids).h], 'descend');
            ids = ids(perm2);
            stepVals = [results(ids).h]';
            errL1 = zeros(numel(ids), 1);
            errLinf = zeros(numel(ids), 1);
            for j = 1:numel(ids)
                R = results(ids(j));
                nRefInterp = interp1(Rref.x, Rref.final_n, R.x, 'pchip');
                nRefInterp = nRefInterp(:);
                diffj = R.final_n - nRefInterp;
                errL1(j) = R.h * sum(abs(diffj));
                errLinf(j) = max(abs(diffj));
            end
            Nref = Rref.N;
            href = Rref.h;
            modeText = sprintf('reference N=%d', Nref);
        end

        orderL1 = compute_orders(stepVals, errL1);
        orderLinf = compute_orders(stepVals, errLinf);

        for j = 1:numel(ids)
            R = results(ids(j));
            gamma_col(end+1,1) = R.gamma; %#ok<AGROW>
            C_col(end+1,1) = R.C; %#ok<AGROW>
            N_col(end+1,1) = R.N; %#ok<AGROW>
            h_col(end+1,1) = R.h; %#ok<AGROW>
            tau_col(end+1,1) = R.tau; %#ok<AGROW>
            N_ref_col(end+1,1) = Nref; %#ok<AGROW>
            h_ref_col(end+1,1) = href; %#ok<AGROW>
            err_L1_T_col(end+1,1) = errL1(j); %#ok<AGROW>
            order_L1_T_col(end+1,1) = orderL1(j); %#ok<AGROW>
            err_Linf_T_col(end+1,1) = errLinf(j); %#ok<AGROW>
            order_Linf_T_col(end+1,1) = orderLinf(j); %#ok<AGROW>
            error_mode_col{end+1,1} = modeText; %#ok<AGROW>
            G_label_col{end+1,1} = R.G_label; %#ok<AGROW>
        end
    end

    if isempty(gamma_col)
        spatialSummary = table();
        return;
    end

    spatialSummary = table(G_label_col, gamma_col, C_col, N_col, h_col, ...
        tau_col, N_ref_col, h_ref_col, err_L1_T_col, order_L1_T_col, ...
        err_Linf_T_col, order_Linf_T_col, error_mode_col, ...
        'VariableNames', {'G','gamma','C','N','h','tau','N_ref','h_ref', ...
        'err_L1_T','order_L1_T','err_Linf_T','order_Linf_T','error_mode'});
end

function caseTable = build_case_table(results)

    if isempty(results)
        caseTable = table();
        return;
    end

    nCases = numel(results);

    G_col = cell(nCases, 1);
    success_col = false(nCases, 1);
    gamma_col = NaN(nCases, 1);
    C_col = NaN(nCases, 1);
    N_col = NaN(nCases, 1);
    h_col = NaN(nCases, 1);
    tau_col = NaN(nCases, 1);
    Nt_col = NaN(nCases, 1);
    err_L1_exact_T_col = NaN(nCases, 1);
    err_Linf_exact_T_col = NaN(nCases, 1);
    rel_mass_change_col = NaN(nCases, 1);
    energy_ratio_col = NaN(nCases, 1);
    min_n_col = NaN(nCases, 1);
    max_n_col = NaN(nCases, 1);
    min_p_col = NaN(nCases, 1);
    max_p_col = NaN(nCases, 1);
    newton_iter_mean_col = NaN(nCases, 1);
    newton_iter_max_col = NaN(nCases, 1);
    newton_res_max_col = NaN(nCases, 1);
    runtime_col = NaN(nCases, 1);
    msg_col = cell(nCases, 1);

    for i = 1:nCases
        R = results(i);
        G_col{i} = R.G_label;
        success_col(i) = R.success;
        gamma_col(i) = R.gamma;
        C_col(i) = R.C;
        N_col(i) = R.N;
        h_col(i) = R.h;
        tau_col(i) = R.tau;
        Nt_col(i) = R.Nt;
        err_L1_exact_T_col(i) = R.err_L1_exact_T;
        err_Linf_exact_T_col(i) = R.err_Linf_exact_T;
        rel_mass_change_col(i) = R.rel_mass_change;
        energy_ratio_col(i) = R.energy_ratio;
        if ~isempty(R.min_n)
            min_n_col(i) = min(R.min_n);
            max_n_col(i) = max(R.max_n);
        end
        if ~isempty(R.min_p)
            min_p_col(i) = min(R.min_p);
            max_p_col(i) = max(R.max_p);
        end
        if ~isempty(R.newton_iter)
            finiteIter = R.newton_iter(isfinite(R.newton_iter));
            if ~isempty(finiteIter)
                newton_iter_mean_col(i) = mean(finiteIter);
                newton_iter_max_col(i) = max(finiteIter);
            end
        end
        if ~isempty(R.newton_res)
            finiteRes = R.newton_res(isfinite(R.newton_res));
            if ~isempty(finiteRes)
                newton_res_max_col(i) = max(finiteRes);
            end
        end
        runtime_col(i) = R.runtime;
        msg_col{i} = R.error_message;
    end

    caseTable = table(G_col, success_col, gamma_col, C_col, N_col, h_col, ...
        tau_col, Nt_col, err_L1_exact_T_col, err_Linf_exact_T_col, ...
        rel_mass_change_col, energy_ratio_col, min_n_col, max_n_col, ...
        min_p_col, max_p_col, newton_iter_mean_col, newton_iter_max_col, ...
        newton_res_max_col, runtime_col, msg_col, ...
        'VariableNames', {'G','success','gamma','C','N','h','tau','Nt', ...
        'err_L1_exact_T','err_Linf_exact_T','rel_mass_change','energy_ratio', ...
        'min_n','max_n','min_p','max_p','Newton_iter_mean','Newton_iter_max', ...
        'Newton_res_max','runtime_sec','error_message'});
end

function idx = find_success_indices(results)
    idx = [];
    for i = 1:numel(results)
        if results(i).success
            idx(end+1) = i; %#ok<AGROW>
        end
    end
end

function orders = compute_orders(step, err)
    orders = NaN(size(err));
    for j = 2:numel(err)
        if err(j) > 0 && err(j-1) > 0 && step(j) > 0 && step(j-1) > 0
            orders(j) = log(err(j-1) / err(j)) / log(step(j-1) / step(j));
        end
    end
end

%% ========================================================================
%% Initial condition, growth functions, and diagnostics
%% ========================================================================

function n = barenblatt_exact(x, t, gamma, C, t0)
% Time-shifted Barenblatt profile from the paper.
%   kappa = gamma/(gamma+1), beta = 1/(gamma+2), s = kappa*(t+t0),
%   n_ex = s^{-beta} * ( C - beta*gamma*x^2/(2*(gamma+1)*s^{2*beta}) )_+^{1/gamma}.

    kappa = gamma / (gamma + 1.0);
    beta = 1.0 / (gamma + 2.0);
    s = kappa * (t + t0);
    if s <= 0
        error('Barenblatt parameter s must be positive. Check t and t0.');
    end
    inside = C - (beta * gamma / (2.0 * (gamma + 1.0))) * (x(:).^2) / (s^(2.0 * beta));
    inside = max(inside, 0.0);
    n = s^(-beta) * inside.^(1.0 / gamma);
end

function C = get_barenblatt_C(params, gamma)
    if abs(gamma - 3) <= 1e-12
        C = params.CForGamma3;
    elseif abs(gamma - 20) <= 1e-12
        C = params.CForGamma20;
    else
        error('This driver is configured for gamma=3 or gamma=20. Received gamma=%g.', gamma);
    end
end

function [Gfun, pH, GpH, GTag, GLabel, exactAvailable, energyPH] = get_growth_case(GCase)

    key = lower(strtrim(GCase));

    switch key
        case {'zero', 'g0', '0', 'none'}
            Gfun = @(p) zeros(size(p));
            pH = 0.0;
            GpH = 0.0;
            GTag = 'G0';
            GLabel = 'G(p)=0';
            exactAvailable = true;
            energyPH = 0.0;

        case {'linear', '1-p', 'one_minus_p', 'g1mp'}
            Gfun = @(p) 1.0 - p;
            pH = 1.0;
            GpH = -1.0;
            GTag = 'G1mp';
            GLabel = 'G(p)=1-p';
            exactAvailable = false;
            energyPH = pH;

        otherwise
            error('Unknown GCase "%s". Use ''zero'' or ''linear''.', GCase);
    end
end

function E = discrete_energy(n, h, gamma, energyPH)
    E = h * sum(n.^(gamma + 1.0) / (gamma + 1.0) - energyPH * n);
end

%% ========================================================================
%% Result structures
%% ========================================================================

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
    R.C = NaN;
    R.t0 = NaN;
    R.G_case = '';
    R.G_tag = '';
    R.G_label = '';
    R.pH = NaN;
    R.GpH = NaN;
    R.energy_pH = NaN;
    R.exact_available = false;
    R.initial_type = '';
    R.output_steps = [];
    R.output_times = [];
    R.snap_n = [];
    R.snap_p = [];
    R.snap_n_exact = [];
    R.snap_files = [];
    R.snap_mass = [];
    R.snap_energy = [];
    R.snap_L1_exact = [];
    R.snap_Linf_exact = [];
    R.final_n = [];
    R.final_p = [];
    R.final_n_exact = [];
    R.time = [];
    R.mass = [];
    R.energy = [];
    R.min_n = [];
    R.max_n = [];
    R.min_p = [];
    R.max_p = [];
    R.rel_mass_change = NaN;
    R.energy_ratio = NaN;
    R.err_L1_exact_T = NaN;
    R.err_Linf_exact_T = NaN;
    R.err_L1_exact_max = NaN;
    R.err_Linf_exact_max = NaN;
    R.newton_iter = [];
    R.newton_res = [];
    R.newton_converged = [];
    R.runtime = NaN;
end

function R = make_failed_result(params, gamma, C, GCase, GTag, GLabel, ...
    pH, energyPH, exactAvailable, NTarget, tauTarget, msg)

    R = make_empty_result();
    R.success = false;
    R.error_message = msg;
    R.a = params.a;
    R.b = params.b;
    R.N = round(NTarget);
    R.h = (params.b - params.a) / round(NTarget);
    R.tau_target = tauTarget;
    R.T = params.T;
    R.gamma = gamma;
    R.C = C;
    R.t0 = params.t0;
    R.G_case = GCase;
    R.G_tag = GTag;
    R.G_label = GLabel;
    R.pH = pH;
    R.energy_pH = energyPH;
    R.exact_available = exactAvailable;
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
%% Plotting
%% ========================================================================

function livePlot = initialize_live_n_plot(x, n, gamma, C, GLabel, N, h, t)

    livePlot.fig = figure('Name', sprintf('Live n, gamma=%g, N=%d', gamma, N));
    set(livePlot.fig, 'Color', 'w');

    livePlot.ax = axes('Parent', livePlot.fig);
    livePlot.line = plot(livePlot.ax, x, n, 'LineWidth', 2.0);

    xlabel(livePlot.ax, '$x$', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel(livePlot.ax, '$n(x,t)$', 'Interpreter', 'latex', 'FontSize', 16);
    title(livePlot.ax, sprintf('%s, $\\gamma=%g$, $C=%g$, $N=%d$, $h=%.3g$, $t=%.5f$', ...
        GLabel, gamma, C, N, h, t), 'Interpreter', 'latex', 'FontSize', 16);

    set(livePlot.ax, 'FontSize', 14, 'LineWidth', 1.0, ...
        'TickLabelInterpreter', 'latex');
    grid(livePlot.ax, 'on');
    box(livePlot.ax, 'on');
    xlim(livePlot.ax, [x(1), x(end)]);
    ylim(livePlot.ax, [0, max(1.0e-12, 1.05 * max(n))]);

    livePlot.gamma = gamma;
    livePlot.C = C;
    livePlot.GLabel = GLabel;
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
        '%s, $\\gamma=%g$, $C=%g$, $N=%d$, $h=%.3g$, $t=%.5f$, ', ...
        'step=%d/%d, Newton=%g, res=%.2e'], ...
        livePlot.GLabel, livePlot.gamma, livePlot.C, livePlot.N, livePlot.h, t, ...
        step, Nt, newtonIter, newtonRes), ...
        'Interpreter', 'latex', 'FontSize', 16);

    drawnow limitrate;
end

function make_case_profile_plot(R, caseDir)
    if ~R.success
        return;
    end

    fig = figure('Visible', 'off', 'Color', 'w');
    hold on; box on;
    plot(R.x, R.final_n, 'o', 'MarkerSize', 4, 'DisplayName', 'Numerical');
    if R.exact_available
        plot(R.x, R.final_n_exact, '-', 'LineWidth', 1.5, 'DisplayName', 'Barenblatt exact');
    end
    xlabel('$x$', 'Interpreter', 'latex');
    ylabel('$n$', 'Interpreter', 'latex');
    title(sprintf('%s, gamma=%g, C=%g, N=%d, tau=%.3g, T=%.3g', ...
        R.G_label, R.gamma, R.C, R.N, R.tau, R.T), 'Interpreter', 'none');
    legend('Location', 'best');
    grid on;
    set(gca, 'FontSize', 12);

    figFile = fullfile(caseDir, sprintf('final_profile_%s_gamma_%s_C_%s_N_%d_tau_%s.png', ...
        R.G_tag, num_tag(R.gamma), num_tag(R.C), R.N, num_tag(R.tau)));
    print(fig, figFile, '-dpng', '-r300');
    close(fig);
end

function make_group_profile_plot(results, groupDir, exactAvailable)
    valid = find_success_indices(results);
    if isempty(valid)
        return;
    end

    fig = figure('Visible', 'off', 'Color', 'w');
    hold on; box on;
    for k = 1:numel(valid)
        R = results(valid(k));
        plot(R.x, R.final_n, 'DisplayName', sprintf('N=%d, tau=%.2g', R.N, R.tau));
    end
    if exactAvailable
        R0 = results(valid(end));
        plot(R0.x, R0.final_n_exact, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Exact');
    end
    xlabel('$x$', 'Interpreter', 'latex');
    ylabel('$n$', 'Interpreter', 'latex');
    title('Final-time profiles', 'Interpreter', 'latex');
    legend('Location', 'best');
    grid on;
    set(gca, 'FontSize', 12);
    print(fig, fullfile(groupDir, 'final_profiles_all_cases.png'), '-dpng', '-r300');
    close(fig);
end

function make_error_summary_plots(temporalSummary, spatialSummary, groupDir, exactAvailable)
    if ~isempty(temporalSummary) && height(temporalSummary) > 0
        fig = figure('Visible', 'off', 'Color', 'w');
        loglog(temporalSummary.tau, temporalSummary.err_L1_T, 'o-', 'LineWidth', 1.5);
        xlabel('$\tau$', 'Interpreter', 'latex');
        ylabel('$L^1$ error at $T$', 'Interpreter', 'latex');
        if exactAvailable
            title('Temporal error against exact Barenblatt', 'Interpreter', 'latex');
        else
            title('Temporal error against numerical reference', 'Interpreter', 'latex');
        end
        grid on; box on;
        print(fig, fullfile(groupDir, 'temporal_error_L1.png'), '-dpng', '-r300');
        close(fig);
    end

    if ~isempty(spatialSummary) && height(spatialSummary) > 0
        fig = figure('Visible', 'off', 'Color', 'w');
        loglog(spatialSummary.h, spatialSummary.err_L1_T, 'o-', 'LineWidth', 1.5);
        xlabel('$h$', 'Interpreter', 'latex');
        ylabel('$L^1$ error at $T$', 'Interpreter', 'latex');
        if exactAvailable
            title('Spatial error against exact Barenblatt', 'Interpreter', 'latex');
        else
            title('Spatial error against numerical reference', 'Interpreter', 'latex');
        end
        grid on; box on;
        print(fig, fullfile(groupDir, 'spatial_error_L1.png'), '-dpng', '-r300');
        close(fig);
    end
end

%% ========================================================================
%% Utilities
%% ========================================================================

function validate_user_parameters(params)
    if params.b <= params.a
        error('Require params.b > params.a.');
    end
    if params.T <= 0
        error('Require params.T > 0.');
    end
    if params.t0 <= 0
        error('Require params.t0 > 0 for the time-shifted Barenblatt profile.');
    end
    if isempty(params.gammaList)
        error('params.gammaList must not be empty.');
    end
    for k = 1:numel(params.gammaList)
        gamma = params.gammaList(k);
        if ~(abs(gamma - 3) <= 1e-12 || abs(gamma - 20) <= 1e-12)
            error('params.gammaList can only contain 3 and/or 20. Received %g.', gamma);
        end
    end
    if isempty(params.NList) || any(params.NList < 2) || any(abs(params.NList - round(params.NList)) > 0)
        error('params.NList must contain positive integers larger than 1.');
    end
    if isempty(params.tauList) || any(params.tauList <= 0)
        error('params.tauList must contain positive time steps.');
    end
    if isempty(params.GCaseList)
        error('params.GCaseList must contain at least one case: ''zero'' or ''linear''.');
    end
    for k = 1:numel(params.GCaseList)
        get_growth_case(params.GCaseList{k});
    end
end

function outputTimes = sanitize_output_times(outputTimes, T)
    outputTimes = outputTimes(:);
    outputTimes = outputTimes(isfinite(outputTimes));
    outputTimes = outputTimes(outputTimes >= 0 & outputTimes <= T);
    outputTimes = [0; outputTimes; T];
    outputTimes = unique(outputTimes, 'stable');
    outputTimes = sort(outputTimes);
end

function list = normalize_cellstr(x)
    if ischar(x)
        list = {x};
    elseif iscell(x)
        list = x(:)';
    else
        error('Expected a character vector or a cell array of character vectors.');
    end
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

function y = round_sig(x, n)
    if nargin < 2
        n = 14;
    end

    y = zeros(size(x));
    idx = (x ~= 0) & isfinite(x);
    if any(idx(:))
        scale = 10.^(n - ceil(log10(abs(x(idx)))));
        y(idx) = round(x(idx) .* scale) ./ scale;
    end
    y(~idx) = x(~idx);
end

function out = strjoin_compat(list, delimiter)
    if isempty(list)
        out = '';
        return;
    end
    out = list{1};
    for k = 2:numel(list)
        out = [out, delimiter, list{k}]; %#ok<AGROW>
    end
end
