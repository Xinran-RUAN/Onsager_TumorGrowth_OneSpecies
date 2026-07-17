clear; clc; close all;

outdir = 'example1_t02';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

% Domain and mesh
a = -20;
b = 20;
h = 1/64;
N = round((b-a)/h);
h = (b-a)/N;
x = a + ((1:N)' - 0.5) * h;

% Time step
tau = 0.01;
T = 100;
Nt = round(T / tau);
time = (0:Nt)' * tau;

% No-growth case
pH = 1.0;
Gfun = @(p) 1 - p;
GpH = -1;

% Two gamma values
gammaList = [3, 20];
CbarList = [1.0, 0.1];
% Parameters in delayed Barenblatt profile
t0 = 2;

opts.tol = 1e-10;
opts.maxit = 50;
opts.verbose = false;

for ig = 1:length(gammaList)

    gamma = gammaList(ig);
    Cbar = CbarList(ig);
    % Initial data from exact Barenblatt profile
    n0 = barenblatt_exact(x, 0, gamma, t0, Cbar);
    n = n0;

    sol_n = zeros(N, Nt+1);
    sol_p = zeros(N, Nt+1);
    err_L1 = zeros(Nt+1, 1);
    mass = zeros(Nt+1, 1);
    energy = zeros(Nt+1, 1);
    newton_iter = zeros(Nt, 1);
    newton_res = zeros(Nt, 1);

    sol_n(:,1) = n;
    sol_p(:,1) = n.^gamma;

    nex = barenblatt_exact(x, 0, gamma, t0, Cbar);
    err_L1(1) = h * sum(abs(n - nex));
    mass(1) = h * sum(n);
    energy(1) = h * sum(n.^(gamma+1) / (gamma+1));

    for k = 1:Nt

        [n, p, F, info] = tumor_step_newton_n(n, h, tau, gamma, pH, Gfun, GpH, opts);

        if ~info.converged
            fprintf('gamma = %d, step %d: Newton not fully converged, residual = %.4e\n', ...
                gamma, k, info.residual);
        end

        sol_n(:,k+1) = n;
        sol_p(:,k+1) = p;

        nex = barenblatt_exact(x, k*tau, gamma, t0, Cbar);
        err_L1(k+1) = h * sum(abs(n - nex));
        mass(k+1) = h * sum(n);
        energy(k+1) = h * sum(n.^(gamma+1) / (gamma+1));

        newton_iter(k) = info.iter;
        newton_res(k) = info.residual;
    end

    datafile = fullfile(outdir, sprintf('example1_gamma%d.mat', gamma));
    save(datafile, ...
        'a', 'b', 'N', 'h', 'x', 'tau', 'T', 'Nt', 'time', ...
        'gamma', 'pH', 't0', 'Cbar', ...
        'sol_n', 'sol_p', 'err_L1', 'mass', 'energy', ...
        'newton_iter', 'newton_res');

    fprintf('Saved %s\n', datafile);
end