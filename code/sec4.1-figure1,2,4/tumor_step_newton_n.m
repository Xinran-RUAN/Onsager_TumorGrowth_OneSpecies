function [nnew, pnew, Fnew, info] = tumor_step_newton_n(nold, h, tau, gamma, pH, Gfun, GpH, opts)
%TUMOR_STEP_NEWTON_N One time step for the fully discrete tumor-growth scheme.
%
% Unknown: density n^{k+1}; pressure is p^{k+1} = (n^{k+1})^gamma.
%
% Fully discrete scheme:
%   (n^{k+1}_i - n^k_i)/tau
%       = (d_h F^{k+1})_i - M_i^k (p_i^{k+1} - pH),
%
%   F_{i+1/2}^{k+1}
%       = n_{i+1/2}^k (p_{i+1}^{k+1} - p_i^{k+1})/h,
%
% with F_{1/2}=F_{N+1/2}=0 and
%   M_i^k = -n_i^k g(p_i^k),
%   g(p)=G(p)/(p-pH), g(pH)=G'(pH).
%
% The nonlinear residual is the unscaled residual used in Appendix A:
%   R_i(u) = u_i - n_i^k + tau*M_i^k*(u_i^gamma-pH)
%            - tau*(d_h F(u))_i.
%
% A damped Newton iteration is used. The line search preserves
% nonnegativity of the Newton trial state.

if nargin < 8 || isempty(opts)
    opts = struct();
end

if ~isfield(opts, 'tol')
    opts.tol = 1e-10;
end
if ~isfield(opts, 'maxit')
    opts.maxit = 1000;
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

% Basic input checks.  The nonnegativity check must be done before forming
% pold=nold.^gamma, otherwise a small negative nold with noninteger gamma
% may create complex values.
if N == 0
    error('nold must be a nonempty vector.');
end
if ~isscalar(h) || ~isfinite(h) || h <= 0
    error('h must be a positive finite scalar.');
end
if ~isscalar(tau) || ~isfinite(tau) || tau <= 0
    error('tau must be a positive finite scalar.');
end
if ~isscalar(gamma) || ~isfinite(gamma) || gamma <= 1
    error('gamma must be a finite scalar larger than 1.');
end
if ~isscalar(pH) || ~isfinite(pH)
    error('pH must be a finite scalar.');
end
if ~isscalar(GpH) || ~isfinite(GpH)
    error('GpH must be a finite scalar.');
end
if any(~isfinite(nold))
    error('nold contains NaN or Inf values.');
end
if min(nold) < 0
    error('The input nold contains negative values.');
end

% Old pressure.
pold = nold.^gamma;

% Interface coefficient A_{i+1/2}^k = n_{i+1/2}^k.
% MATLAB indexing:
%   A(1)   = A_{1/2}
%   A(j)   = A_{j-1/2}
%   A(N+1) = A_{N+1/2}
A = zeros(N+1, 1);
A(2:N) = 0.5 * (nold(1:N-1) + nold(2:N));
A(1) = 0.0;
A(N+1) = 0.0;

% Reaction mobility M_i^k = -n_i^k g(p_i^k).
M = zeros(N, 1);
idx = abs(pold - pH) > opts.tolH * (1 + abs(pH));

if any(idx)
    Gvals = Gfun(pold(idx));
    Gvals = Gvals(:);
    if numel(Gvals) ~= nnz(idx)
        error('Gfun must return an array with the same number of entries as its input.');
    end
    if any(~isfinite(Gvals))
        error('Gfun returned NaN or Inf values.');
    end
    M(idx) = -nold(idx) .* Gvals ./ (pold(idx) - pH);
end

% Limiting value when pold is close to pH:
%   M = - n^k G'(pH)
M(~idx) = -nold(~idx) .* GpH;

% Initial guess u^(0)=n^k.
n = nold;

R = compute_residual(n, nold, A, M, h, tau, gamma, pH);
Rinf = norm(R, inf);

info.converged = false;
info.iter = 0;
info.residual = Rinf;
info.residual2 = norm(R, 2);
info.lambda = [];

if opts.verbose
    fprintf('Newton iteration: initial residual_inf = %.4e\n', Rinf);
end

for it = 1:opts.maxit

    R = compute_residual(n, nold, A, M, h, tau, gamma, pH);
    resInf = norm(R, inf);

    if resInf <= opts.tol * (1 + norm(nold, inf))
        info.converged = true;
        info.iter = it - 1;
        info.residual = resInf;
        info.residual2 = norm(R, 2);
        break;
    end

    J = compute_jacobian(n, A, M, h, tau, gamma);

    dn = -(J \ R);

    % Backtracking line search.  The Armijo condition uses the 2-norm as in
    % Appendix A, while the stopping criterion above uses the infinity norm.
    lambda = 1.0;
    accepted = false;
    res2 = norm(R, 2);

    for ls = 1:opts.lineSearchMaxit

        ntrial = n + lambda * dn;

        % Preserve nonnegativity of the Newton trial state.
        if min(ntrial) >= 0 && all(isfinite(ntrial))

            Rtrial = compute_residual(ntrial, nold, A, M, h, tau, gamma, pH);
            trialNorm2 = norm(Rtrial, 2);

            if isfinite(trialNorm2) && trialNorm2 <= (1 - opts.armijo * lambda) * res2
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
        info.residual = resInf;
        info.residual2 = res2;
        break;
    end

    n = ntrial;

    Raccepted = Rtrial;
    info.lambda = [info.lambda; lambda]; %#ok<AGROW>
    info.iter = it;
    info.residual = norm(Raccepted, inf);
    info.residual2 = trialNorm2;

    if opts.verbose
        fprintf('  it = %2d, residual_inf = %.4e, lambda = %.3e\n', ...
                it, info.residual, lambda);
    end
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
    info.residual2 = norm(R, 2);
end

end


function R = compute_residual(n, nold, A, M, h, tau, gamma, pH)
% Unscaled residual from Appendix A:
%   R_i = n_i - nold_i + tau*M_i*(p_i - pH) - tau*(d_h F)_i.

N = length(n);
p = n.^gamma;

F = zeros(N+1, 1);
F(2:N) = A(2:N) .* (p(2:N) - p(1:N-1)) / h;
F(1) = 0.0;
F(N+1) = 0.0;

dF = (F(2:N+1) - F(1:N)) / h;

R = n - nold + tau * M .* (p - pH) - tau * dF;

end


function J = compute_jacobian(n, A, M, h, tau, gamma)
% Jacobian of the unscaled residual with respect to n.

N = length(n);
ng = gamma * n.^(gamma - 1);

main = 1.0 ...
     + tau * (A(1:N) + A(2:N+1)) .* ng / h^2 ...
     + tau * M .* ng;

lower = -tau * A(2:N) .* ng(1:N-1) / h^2;
upper = -tau * A(2:N) .* ng(2:N) / h^2;

J = spdiags(main, 0, N, N) ...
  + sparse(2:N, 1:N-1, lower, N, N) ...
  + sparse(1:N-1, 2:N, upper, N, N);

end
