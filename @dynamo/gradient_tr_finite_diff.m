function ret = gradient_tr_finite_diff(self, t, k, c)
% Gradient of error_tr by finite difference method.

% f'(x) = (f(x + eps) -f(x))/eps
% Trivial and relatively slow, but a good reference point.
%
% Uses H{t}, U{t}, L{t+1} and g


[P, epsilon] = self.finite_diff_P(t, k, c);

% We may compute the finite diff approximation at three locations:
if 0
    % at P
    dPdu = (P -self.cache.P{t, k}) / epsilon;
    dgdu = partial_trace(self.cache.L{t+1, k} * (dPdu * self.cache.U{t, k}), self.system.dimSE, 1);
else
    g = partial_trace(self.cache.L{t+1, k} * (P * self.cache.U{t, k}), self.system.dimSE, 1);
    if 0
        % at g
        dgdu = (g -self.cache.g{k}) / epsilon;
    else
        % at E
        E = self.config.f_max -trace(svd(g));
        ret = (E -self.cache.E) / epsilon;
        return
    end
end
ret = -trace_matmul(self.cache.VUdagger, dgdu);
