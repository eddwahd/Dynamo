function ret = gradient_tr_finite_diff(self, t, k, c)
% Gradient of error_tr by finite difference method.

% E'(x) = (E(x + eps) - E(x))/eps
% Trivial and relatively slow, but a good reference point.
%
% Uses H{t}, U{t}, L{t+1} and g

[P, epsilon] = self.finite_diff_P(t, k, c);
g_at_eps = partial_trace(self.cache.L{t+1, k} * (P * self.cache.U{t, k}), self.system.dimSE, 1);

% just do dQdr by finite diff
dQdr = (g_at_eps -self.cache.g{k}) / epsilon;
ret = -trace_matmul(self.cache.VUdagger, dQdr);

% TODO alternative: finite diff all the way
%[U, S, V] = svd(g_at_eps);
%E_at_eps = self.config.f_max -trace(S);
%ret = (E_at_eps -self.cache.E) / epsilon;
