% Parameters
a = -5;
b = 5;
N = 1000;
h = (b - a) / N;
x = a + ((1:N)' - 0.5) * h;

tau = 1e-2;
T = 2.;
Nt = round(T / tau);

gamma = 100;
pH = 1.0;
alpha = 1.0;

Gfun = @(p) alpha * (pH - p);
GpH = -alpha;

% Initial data: compactly supported patch
R0 = 1.0;
n0 = 0.8 * double(abs(x) <= R0);
n = n0;

opts.tol = 1e-10;
opts.maxit = 50;
opts.verbose = false;

% Storage
sol_n = zeros(N, Nt+1);
sol_p = zeros(N, Nt+1);
sol_n(:,1) = n;
sol_p(:,1) = n.^gamma;

% Plot setting
plotEvery = round(0.02 / tau);
plotEvery = max(plotEvery, 1);

figure;
subplot(2,1,1);
hN = plot(x, n, 'LineWidth', 2);
xlabel('x');
ylabel('n');
title(sprintf('Density, t = %.4f', 0));
ylim([-0.05, 1.1]);
grid on;

subplot(2,1,2);
hP = plot(x, n.^gamma, 'LineWidth', 2);
xlabel('x');
ylabel('p');
title(sprintf('Pressure, t = %.4f', 0));
ylim([-0.05, 1.1]);
grid on;

drawnow;

for k = 1:Nt

    [n, p, F, info] = tumor_step_newton_n(n, h, tau, gamma, pH, Gfun, GpH, opts);

    if ~info.converged
        fprintf('Step %d: Newton not fully converged, residual = %.4e\n', ...
                k, info.residual);
    end

    sol_n(:,k+1) = n;
    sol_p(:,k+1) = p;

    if mod(k, plotEvery) == 0 || k == Nt
        subplot(2,1,1);
        set(hN, 'YData', n);
        title(sprintf('Density, t = %.4f', k * tau));

        subplot(2,1,2);
        set(hP, 'YData', p);
        title(sprintf('Pressure, t = %.4f', k * tau));

        drawnow;
    end
end