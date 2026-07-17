% plot_barenblatt_initial_t0.m
% Compare Barenblatt initial profiles n0(x)=n_ex(x,0)
% for different choices of t0.

clear; clc; close all;

%% Parameters
a = -5;
b =  5;
Nx = 4000;
x = linspace(a, b, Nx);

% Choose gamma = 3 or gamma = 20
gamma = 20;
% gamma = 20;

% Choose C according to the paper
if gamma == 3
    C = 1.0;
elseif gamma == 20
    C = 0.1;
else
    C = 1.0;   % default, can be changed manually
end

% Different time shifts to compare
t0List = [1e-2, 5e-2, 1e-1, 1, 2, 3];

% Barenblatt parameters
kappa = gamma / (gamma + 1);
beta  = 1 / (gamma + 2);

%% Storage
nMat = zeros(length(t0List), Nx);
pMat = zeros(length(t0List), Nx);
R0List = zeros(length(t0List), 1);
maxNList = zeros(length(t0List), 1);
maxPList = zeros(length(t0List), 1);
massList = zeros(length(t0List), 1);

%% Compute initial profiles
for m = 1:length(t0List)
    t0 = t0List(m);

    s = kappa * t0;

    inside = C - beta * gamma / (2 * (gamma + 1)) * x.^2 ./ (s.^(2*beta));
    inside = max(inside, 0);

    n0 = s^(-beta) * inside.^(1/gamma);
    p0 = n0.^gamma;

    nMat(m,:) = n0;
    pMat(m,:) = p0;

    % Theoretical support radius
    R0List(m) = sqrt(2 * (gamma + 1) * C / (beta * gamma)) * s^beta;

    maxNList(m) = max(n0);
    maxPList(m) = max(p0);
    massList(m) = trapz(x, n0);
end

%% Print diagnostics
fprintf('gamma = %g, C = %g\n', gamma, C);
fprintf('---------------------------------------------------------------\n');
fprintf('   t0          support R0       max(n0)        max(p0)       mass\n');
fprintf('---------------------------------------------------------------\n');
for m = 1:length(t0List)
    fprintf('%10.3e   %12.5e   %12.5e   %12.5e   %12.5e\n', ...
        t0List(m), R0List(m), maxNList(m), maxPList(m), massList(m));
end
fprintf('---------------------------------------------------------------\n');

%% Plot density profiles
figure;
hold on; box on; grid on;

for m = 1:length(t0List)
    plot(x, nMat(m,:), 'LineWidth', 1.5);
end

xlabel('x');
ylabel('n_0(x)');
title(sprintf('Barenblatt initial density, \\gamma = %g, C = %g', gamma, C));

legText = cell(length(t0List), 1);
for m = 1:length(t0List)
    legText{m} = sprintf('t_0 = %.0e', t0List(m));
end
legend(legText, 'Location', 'best');
set(gca, 'FontSize', 12);

%% Plot pressure profiles
figure;
hold on; box on; grid on;

for m = 1:length(t0List)
    plot(x, pMat(m,:), 'LineWidth', 1.5);
end

xlabel('x');
ylabel('p_0(x)=n_0(x)^\gamma');
title(sprintf('Initial pressure, \\gamma = %g, C = %g', gamma, C));

legend(legText, 'Location', 'best');
set(gca, 'FontSize', 12);

%% Plot zoomed density near the interface
figure;
hold on; box on; grid on;

for m = 1:length(t0List)
    plot(x, nMat(m,:), 'LineWidth', 1.5);
end

xlabel('x');
ylabel('n_0(x)');
title('Zoomed density profile');
legend(legText, 'Location', 'best');
set(gca, 'FontSize', 12);

xlim([-6, 6]);