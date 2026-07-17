function [nnew, pnew, Fnew, info] = tumor_step_newton_n(nold, h, tau, gamma, pH, Gfun, GpH, opts)
%TUMOR_STEP_NEWTON_N One time step for the fully discrete tumor-growth scheme.
%
% The main unknown is n^{k+1}. The pressure is p^{k+1} = (n^{k+1})^gamma.
%
% Scheme:
%   (n^{k+1}_i - n^k_i)/tau
%     = (d_h F^{k+1})_i - M_i^k (p_i^{k+1} - pH),
%
%   F_{i+1/2}^{k+1}
%     = n_{i+1/2}^k (p_{i+1}^{k+1} - p_i^{k+1})/h,
%
% with F_{1/2}=F_{N+1/2}=0.
%
% Inputs:
%   nold  : density at time step k, column vector of length N
%   h     : mesh size
%   tau   : time step
%   gamma : exponent in p = n^gamma
%   pH    : homeostatic pressure
%   Gfun  : function handle for G(p)
%   GpH   : value G'(pH), used for the limiting definition of M
%   opts  : optional structure
%
% Outputs:
%   nnew  : density at time step k+1
%   pnew  : pressure at time step k+1
%   Fnew  : auxiliary interface quantity, length N+1
%   info  : iteration information

if nargin < 8
    opts = struct();
end

if ~isfield(opts, 'tol')
    opts.tol = 1e-9;
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

nold = nold(:);
N = length(nold);

% Old pressure
pold = nold.^gamma;

% Interface coefficient A_{i+1/2}^k = n_{i+1/2}^k.
% In MATLAB indexing:
%   A(1)   = A_{1/2}
%   A(j)   = A_{j-1/2}
%   A(N+1) = A_{N+1/2}
A = zeros(N+1, 1);
A(2:N) = 0.5 * (nold(1:N-1) + nold(2:N));
A(1) = 0.0;
A(N+1) = 0.0;

% Compute M_i^k.
M = zeros(N, 1);
idx = abs(pold - pH) > opts.tolH * (1 + abs(pH));

M(idx) = -nold(idx) .* Gfun(pold(idx)) ./ (pold(idx) - pH);

% Limiting value when pold is close to pH:
%   M = - n^k G'(pH)
M(~idx) = -nold(~idx) .* GpH;

% Initial guess
n = nold;

% If the initial data contains very small negative values due to roundoff,
% stop rather than modifying the scheme silently.
if min(n) < 0
    error('The input nold contains negative values.');
end

R = compute_residual(n, nold, A, M, h, tau, gamma, pH);
R0 = norm(R, inf);

info.converged = false;
info.iter = 0;
info.residual = R0;
info.lambda = [];

if opts.verbose
    fprintf('Newton iteration: initial residual = %.4e\n', R0);
end

for it = 1:opts.maxit

    R = compute_residual(n, nold, A, M, h, tau, gamma, pH);
    resnorm = norm(R, inf);

    if resnorm <= opts.tol * (1 + norm(nold, inf))
        info.converged = true;
        info.iter = it - 1;
        info.residual = resnorm;
        break;
    end

    J = compute_jacobian(n, A, M, h, tau, gamma);

    dn = -J \ R;

    % Backtracking line search.
    lambda = 1.0;
    accepted = false;

    for ls = 1:opts.lineSearchMaxit

        ntrial = n + lambda * dn;

        % Enforce nonnegativity of the Newton trial state.
        if min(ntrial) >= 0

            Rtrial = compute_residual(ntrial, nold, A, M, h, tau, gamma, pH);
            trialNorm = norm(Rtrial, inf);

            % Armijo-type residual decrease condition.
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

    if opts.verbose
        fprintf('  it = %2d, residual = %.4e, lambda = %.3e\n', ...
                it, trialNorm, lambda);
    end

    info.lambda = [info.lambda; lambda];
    info.iter = it;
    info.residual = trialNorm;
end

nnew = n;
pnew = nnew.^gamma;

Fnew = zeros(N+1, 1);
Fnew(2:N) = A(2:N) .* (pnew(2:N) - pnew(1:N-1)) / h;
Fnew(1) = 0.0;
Fnew(N+1) = 0.0;

if ~info.converged
    R = compute_residual(nnew, nold, A, M, h, tau, gamma, pH);
    info.residual = norm(R, inf);
end

end


function R = compute_residual(n, nold, A, M, h, tau, gamma, pH)
% Residual:
%   R_i = (n_i - nold_i)/tau - (d_h F)_i + M_i (p_i - pH).

N = length(n);
p = n.^gamma;

F = zeros(N+1, 1);
F(2:N) = A(2:N) .* (p(2:N) - p(1:N-1)) / h;
F(1) = 0.0;
F(N+1) = 0.0;

dF = (F(2:N+1) - F(1:N)) / h;

R = (n - nold) / tau - dF + M .* (p - pH);

end


function J = compute_jacobian(n, A, M, h, tau, gamma)
% Jacobian of the residual with respect to n.

N = length(n);

ng = gamma * n.^(gamma - 1);

main = 1/tau ...
     + (A(1:N) + A(2:N+1)) .* ng / h^2 ...
     + M .* ng;

lower = -A(2:N) .* ng(1:N-1) / h^2;
upper = -A(2:N) .* ng(2:N) / h^2;

J = spdiags(main, 0, N, N) ...
  + sparse(2:N, 1:N-1, lower, N, N) ...
  + sparse(1:N-1, 2:N, upper, N, N);

end

