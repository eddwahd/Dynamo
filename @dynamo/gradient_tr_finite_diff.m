function ret = gradient_tr_finite_diff(self, t, k, c)
% Gradient of error_tr by finite difference method.

% E'(x) = (E(x + eps) - E(x))/eps
% Trivial and relatively slow, but a good reference point.
%
% Uses H{t}, U{t}, L{t+1} and g

epsilon = self.config.epsilon;
if c < 0
    % modify tau, re-integrate bin (plus the ones following it!)
    % TODO NOTE with crosstalk and varying taus, all the bins after the one whose
    % length we change will have different propagators because the
    % phases of the slowly rotating terms depend on (absolute) time! keep taus constant?
    tau_eps = self.seq.tau(t) +self.seq.tau_deriv(t) * epsilon;
    P_epsilon = expm(-tau_eps * self.cache.H{t, k});
else
    H_eps = self.cache.H{t, k} +(epsilon * self.seq.fields_deriv(t, c)) * self.system.B{c, k};
    P_epsilon = expm(-self.seq.tau(t) * H_eps);
    %cache.with_modified_control(t,c, f'(t,c)).integrate_bin(t, k);
    % f'(t,c) \approx f(t,c) +epsilon * self.seq.fields_deriv(t, c)
end
g_at_eps_point = partial_trace(self.cache.L{t+1, k} * (P_epsilon * self.cache.U{t, k}), self.system.dimSE, 1);

% just do dQdr by finite diff
dQdr = (g_at_eps_point -self.cache.g{k}) / epsilon;
ret = -trace_matmul(self.cache.VUdagger, dQdr);

% TODO alternative: finite diff all the way
%[U, S, V] = svd(g_at_eps_point);
%E_at_eps_point = self.config.f_max -trace(S);
%ret = (E_at_eps_point -self.cache.E) / epsilon;
