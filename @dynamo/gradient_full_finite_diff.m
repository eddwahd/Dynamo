function ret = gradient_full_finite_diff(self, t, k, c)
% Gradient of error_full by finite difference method.

% E'(x) = (E(x + eps) - E(x))/eps
% Trivial and relatively slow, but a good reference point.
%
% Uses H{t}, U{t} and L{t+1}.

[P, epsilon] = self.finite_diff_P(t, k, c);
X_S = partial_trace(self.cache.L{t+1, k} * (P * self.cache.U{t, k}), self.system.dimSE, 2);

E_at_eps = 0.5 * norm2(self.system.X_final -X_S);
ret = (E_at_eps -self.cache.E) / epsilon;
