function ret = gradient_full_finite_diff(self, t, c)
% Gradient of error_full by finite difference method.

% E'(x) = (E(x + eps) - E(x))/eps
% Trivial and relatively slow, but a good reference point.
%
% Uses H{t}, U{t} and L{t+1}.

epsilon = self.config.epsilon;
if c < 0
    tau_eps = self.seq.tau(t) +self.seq.tau_deriv(t) * epsilon;
    P_epsilon = expm(-tau_eps * self.cache.H{t});
else
    H_eps = self.cache.H{t} +(epsilon * self.seq.fields_deriv(t, c)) * self.system.B{c};
    P_epsilon = expm(-self.seq.tau(t) * H_eps);
end

X_S = partial_trace(self.cache.L{t+1} * (P_epsilon * self.cache.U{t}), self.system.dimSE, 2);

E_at_eps_point = 0.5 * norm2(self.system.X_final -X_S) / self.system.norm2;
ret = (E_at_eps_point - self.cache.E) / epsilon;
