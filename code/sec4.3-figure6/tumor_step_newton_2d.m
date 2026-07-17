function [nnew, pnew, info] = tumor_step_newton_2d(nold, h, tau, gamma, pH, Gfun, GpH, opts)
%TUMOR_STEP_NEWTON_2D One time step of the 2D fully discrete scheme.
%
% Main unknown is n^{k+1}. Pressure is p^{k+1} = (n^{k+1})^gamma.

if nargin < 8
    opts = struct();
end
if ~isfield(opts, 'tol')
    opts.tol = 1e-10;
end
if ~isfield(opts, 'maxit')
    opts.maxit = 50;
end
if ~isfield(opts, 'lineSearchMaxit')
    opts.lineSearchMaxit = 40;
end
if ~isfield(opts, 'tolH')
    opts.tolH = 1e-12;
end
if ~isfield(opts, 'armijo')
    opts.armijo = 1e-4;
end
if ~isfield(opts, 'verbose')
    opts.verbose = false;
end

[Ny, Nx] = size(nold);

if min(nold(:)) < 0
    error('The input nold contains negative values.');
end

pold = nold.^gamma;

% Interface coefficients in x direction.
% Ax has size Ny x (Nx+1).
% Ax(:,1) and Ax(:,Nx+1) are boundary interfaces.
Ax = zeros(Ny, Nx+1);
Ax(:,2:Nx) = 0.5 * (nold(:,1:Nx-1) + nold(:,2:Nx));

% Interface coefficients in y direction.
% Ay has size (Ny+1) x Nx.
% Ay(1,:) and Ay(Ny+1,:) are boundary interfaces.
Ay = zeros(Ny+1, Nx);
Ay(2:Ny,:) = 0.5 * (nold(1:Ny-1,:) + nold(2:Ny,:));

% Compute M^k.
M = zeros(Ny, Nx);
idx = abs(pold - pH) > opts.tolH * (1 + abs(pH));

M(idx) = -nold(idx) .* Gfun(pold(idx)) ./ (pold(idx) - pH);
M(~idx) = -nold(~idx) .* GpH;

% Initial guess
n = nold;

R = compute_residual_2d(n, nold, Ax, Ay, M, h, tau, gamma, pH);
res0 = norm(R(:), inf);

info.converged = false;
info.iter = 0;
info.residual = res0;
info.lambda = [];

for it = 1:opts.maxit

    R = compute_residual_2d(n, nold, Ax, Ay, M, h, tau, gamma, pH);
    resnorm = norm(R(:), inf);

    if resnorm <= opts.tol * (1 + norm(nold(:), inf))
        info.converged = true;
        info.iter = it - 1;
        info.residual = resnorm;
        break;
    end

    J = compute_jacobian_2d(n, Ax, Ay, M, h, tau, gamma);

    dn_vec = -J \ R(:);
    dn = reshape(dn_vec, Ny, Nx);

    lambda = 1.0;
    accepted = false;

    for ls = 1:opts.lineSearchMaxit

        ntrial = n + lambda * dn;

        if min(ntrial(:)) >= 0
            Rtrial = compute_residual_2d(ntrial, nold, Ax, Ay, M, h, tau, gamma, pH);
            trialNorm = norm(Rtrial(:), inf);

            if trialNorm <= (1 - opts.armijo * lambda) * resnorm
                accepted = true;
                break;
            end
        end

        lambda = 0.5 * lambda;
    end

    if ~accepted
        warning('Damped Newton line search failed at iteration %d.', it);
        info.converged = false;
        info.iter = it;
        info.residual = resnorm;
        break;
    end

    n = ntrial;

    info.lambda = [info.lambda; lambda];
    info.iter = it;
    info.residual = trialNorm;

    if opts.verbose
        fprintf('  Newton it = %d, residual = %.4e, lambda = %.3e\n', ...
                it, trialNorm, lambda);
    end
end

nnew = n;
pnew = nnew.^gamma;

if ~info.converged
    R = compute_residual_2d(nnew, nold, Ax, Ay, M, h, tau, gamma, pH);
    info.residual = norm(R(:), inf);
end

end


function R = compute_residual_2d(n, nold, Ax, Ay, M, h, tau, gamma, pH)

[Ny, Nx] = size(n);

p = n.^gamma;

Fx = zeros(Ny, Nx+1);
Fy = zeros(Ny+1, Nx);

Fx(:,2:Nx) = Ax(:,2:Nx) .* (p(:,2:Nx) - p(:,1:Nx-1)) / h;
Fy(2:Ny,:) = Ay(2:Ny,:) .* (p(2:Ny,:) - p(1:Ny-1,:)) / h;

divF = (Fx(:,2:Nx+1) - Fx(:,1:Nx)) / h ...
     + (Fy(2:Ny+1,:) - Fy(1:Ny,:)) / h;

R = (n - nold) / tau - divF + M .* (p - pH);

end


function J = compute_jacobian_2d(n, Ax, Ay, M, h, tau, gamma)

[Ny, Nx] = size(n);
Ntot = Nx * Ny;

D = gamma * n.^(gamma - 1);
ind = reshape(1:Ntot, Ny, Nx);

h2 = h^2;

main = 1/tau ...
     + D .* (Ax(:,1:Nx) + Ax(:,2:Nx+1) ...
            + Ay(1:Ny,:) + Ay(2:Ny+1,:)) / h2 ...
     + M .* D;

rows = ind(:);
cols = ind(:);
vals = main(:);

% Left neighbor
if Nx > 1
    r = ind(:,2:Nx);
    c = ind(:,1:Nx-1);
    v = -Ax(:,2:Nx) .* D(:,1:Nx-1) / h2;

    rows = [rows; r(:)];
    cols = [cols; c(:)];
    vals = [vals; v(:)];

    % Right neighbor
    r = ind(:,1:Nx-1);
    c = ind(:,2:Nx);
    v = -Ax(:,2:Nx) .* D(:,2:Nx) / h2;

    rows = [rows; r(:)];
    cols = [cols; c(:)];
    vals = [vals; v(:)];
end

% Lower neighbor in y
if Ny > 1
    r = ind(2:Ny,:);
    c = ind(1:Ny-1,:);
    v = -Ay(2:Ny,:) .* D(1:Ny-1,:) / h2;

    rows = [rows; r(:)];
    cols = [cols; c(:)];
    vals = [vals; v(:)];

    % Upper neighbor in y
    r = ind(1:Ny-1,:);
    c = ind(2:Ny,:);
    v = -Ay(2:Ny,:) .* D(2:Ny,:) / h2;

    rows = [rows; r(:)];
    cols = [cols; c(:)];
    vals = [vals; v(:)];
end

J = sparse(rows, cols, vals, Ntot, Ntot);

end