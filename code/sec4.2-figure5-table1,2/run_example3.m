clear; clc; close all;

outdir = 'example3';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

% Parameters
a = -5;
b = 5;
N = 6400;
h = (b - a) / N;
x = a + ((1:N)' - 0.5) * h;

tau = 1e-4;
T = 1.0;
Nt = round(T / tau);
time = (0:Nt)' * tau;

pH = 1.0;
alpha = 1.0;

% Growth function G(p)=1-p
Gfun = @(p) alpha * (pH - p);
GpH = -alpha;

% Hele-Shaw initial radius
R0 = 1.0;

% Gamma values
gammaList = [640, 1280];

% Threshold for numerical front location
theta = 0.5;

opts.tol = 1e-8;
opts.maxit = 100;
opts.verbose = false;

for ig = 1:length(gammaList)

    gamma = gammaList(ig);

    fprintf('Running Example 3 with gamma = %d\n', gamma);

    % Well-prepared initial data
    p0_hs = 1 - cosh(x) / cosh(R0);
    p0_hs = max(p0_hs, 0);

    n0 = p0_hs.^(1/gamma);
    n = n0;
    p = n.^gamma;

    % Storage
    sol_n = zeros(N, Nt+1);
    sol_p = zeros(N, Nt+1);

    R_num = zeros(Nt+1, 1);
    R_exact = zeros(Nt+1, 1);

    err_n_L1 = zeros(Nt+1, 1);
    err_p_L1 = zeros(Nt+1, 1);

    max_n = zeros(Nt+1, 1);
    max_p = zeros(Nt+1, 1);
    min_n = zeros(Nt+1, 1);

    energy = zeros(Nt+1, 1);

    newton_iter = zeros(Nt, 1);
    newton_res = zeros(Nt, 1);

    % Initial storage
    sol_n(:,1) = n;
    sol_p(:,1) = p;

    R_exact(1) = asinh(exp(0) * sinh(R0));
    p_hs = 1 - cosh(x) / cosh(R_exact(1));
    p_hs = max(p_hs, 0);
    n_hs = double(abs(x) <= R_exact(1));

    idx_front = find(n > theta);
    if isempty(idx_front)
        R_num(1) = 0;
    else
        R_num(1) = max(abs(x(idx_front)));
    end

    err_n_L1(1) = h * sum(abs(n - n_hs));
    err_p_L1(1) = h * sum(abs(p - p_hs));

    min_n(1) = min(n);
    max_n(1) = max(n);
    max_p(1) = max(p);
    energy(1) = h * sum(n.^(gamma+1) / (gamma+1) - pH * n);

    % Plot during computation
    plotEvery = max(round(0.02 / tau), 1);

    figure(100 + ig);
    subplot(2,1,1);
    hN = plot(x, n, 'LineWidth', 2);
    xlabel('x');
    ylabel('n');
    title(sprintf('Example 3, gamma = %d, density, t = %.4f', gamma, 0));
    ylim([-0.05, 1.1]);
    grid on;

    subplot(2,1,2);
    hP = plot(x, p, 'LineWidth', 2);
    xlabel('x');
    ylabel('p');
    title(sprintf('Example 3, gamma = %d, pressure, t = %.4f', gamma, 0));
    ylim([-0.05, 1.1]);
    grid on;
    drawnow;

    for k = 1:Nt

        [n, p, F, info] = tumor_step_newton_n(n, h, tau, gamma, pH, Gfun, GpH, opts);

        if ~info.converged
            fprintf('gamma = %d, step %d: Newton not fully converged, residual = %.4e\n', ...
                    gamma, k, info.residual);
        end

        sol_n(:,k+1) = n;
        sol_p(:,k+1) = p;

        t = k * tau;

        % Exact Hele-Shaw solution
        R_exact(k+1) = asinh(exp(t) * sinh(R0));
        p_hs = 1 - cosh(x) / cosh(R_exact(k+1));
        p_hs = max(p_hs, 0);
        n_hs = double(abs(x) <= R_exact(k+1));

        % Numerical front location
        idx_front = find(n > theta);
        if isempty(idx_front)
            R_num(k+1) = 0;
        else
            R_num(k+1) = max(abs(x(idx_front)));
        end

        % Errors
        err_n_L1(k+1) = h * sum(abs(n - n_hs));
        err_p_L1(k+1) = h * sum(abs(p - p_hs));

        % Diagnostics
        min_n(k+1) = min(n);
        max_n(k+1) = max(n);
        max_p(k+1) = max(p);
        energy(k+1) = h * sum(n.^(gamma+1) / (gamma+1) - pH * n);

        newton_iter(k) = info.iter;
        newton_res(k) = info.residual;

        if mod(k, plotEvery) == 0 || k == Nt
            subplot(2,1,1);
            set(hN, 'YData', n);
            title(sprintf('Example 3, gamma = %d, density, t = %.4f', gamma, t));

            subplot(2,1,2);
            set(hP, 'YData', p);
            title(sprintf('Example 3, gamma = %d, pressure, t = %.4f', gamma, t));

            drawnow;
        end
    end

    datafile = fullfile(outdir, sprintf('example3_gamma%d_N%d.mat', gamma, N));

    save(datafile, ...
        'a', 'b', 'N', 'h', 'x', 'tau', 'T', 'Nt', 'time', ...
        'gamma', 'pH', 'alpha', 'R0', 'theta', ...
        'sol_n', 'sol_p', ...
        'R_num', 'R_exact', ...
        'err_n_L1', 'err_p_L1', ...
        'min_n', 'max_n', 'max_p', ...
        'energy', 'newton_iter', 'newton_res');

    fprintf('Saved %s\n', datafile);
    fprintf('Final L1 error for n = %.6e\n', err_n_L1(end));
    fprintf('Final L1 error for p = %.6e\n', err_p_L1(end));
end

plot_example3_profiles;
plot_example3_radius;
plot_example3_errors;

fprintf('Example 3 finished.\n');