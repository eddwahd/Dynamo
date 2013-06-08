function ret = gradient_g_finite_diff(self, t, k, c)
% Gradient of the auxiliary function g by finite difference method.
% Special case of gradient_tr_finite_diff.

% f'(x) = (f(x + eps) -f(x))/eps
% Trivial and relatively slow, but a good reference point.
%
% Uses H{t}, U{t}, L{t+1} and g.


[P, epsilon] = self.finite_diff_P(t, k, c);

% We may compute the finite diff approximation at three locations:
if 0
    % at P
    dPdu = (P -self.cache.P{t, k}) / epsilon;
    dgdu = trace_matmul(self.cache.L{t+1, k}, dPdu * self.cache.U{t, k});
else
    g = trace_matmul(self.cache.L{t+1, k}, P * self.cache.U{t, k});
    if 0
        % at g
        dgdu = (g -self.cache.g{k}) / epsilon;
    else
        % at E
        E = self.config.f_max -abs(g);
        ret = (E -self.cache.E) / epsilon;
        return
    end
end
ret = -self.cache.VUdagger * dgdu;
