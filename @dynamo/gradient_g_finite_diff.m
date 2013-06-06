function ret = gradient_g_finite_diff(self, t, k, c)
% Gradient of the auxiliary function g by finite difference method.
% Special case of gradient_tr_finite_diff.

% g'(x) = (g(x + eps) - g(x))/eps
% Trivial and relatively slow, but a good reference point.
%
% Uses g, H{t}, U{t} and L{t+1}.

[P, epsilon] = self.finite_diff_P(t, k, c);
g_at_eps = trace_matmul(self.cache.L{t+1, k}, P * self.cache.U{t, k});

ret = (g_at_eps -self.cache.g{k}) / epsilon;
ret = ret * -self.cache.VUdagger;
