function n = barenblatt_exact(x, t, gamma, t0, Cbar)
%BARENBLATT_EXACT Delayed Barenblatt solution for
%
%   n_t = kappa * (n^(gamma+1))_{xx},
%
% where kappa = gamma/(gamma+1).
%
% This corresponds to
%
%   n_t = (n (n^gamma)_x)_x.
%
% The formula is the standard one-dimensional Barenblatt profile
% with the rescaled time s = kappa * (t + t0).

kappa = gamma / (gamma + 1);
beta = 1 / (gamma + 2);

s = kappa * (t + t0);

coef = beta * gamma / (2 * (gamma + 1));

inside = Cbar - coef * x.^2 ./ (s.^(2 * beta));
inside = max(inside, 0);

n = s.^(-beta) .* inside.^(1/gamma);
end
