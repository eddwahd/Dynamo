function ret = gradient_open_finite_diff(self, t, c, epsilon)
% Gradient of error_open by finite difference method.

% E_open'(x) = (E_open(x + eps) - E_open(x))/eps
% Trivial and relatively slow, but a good reference point.
%
% Uses H{t}, U{t} and L{t+1}.

if c < 0
    tau_eps = self.seq.tau(t) +self.seq.tau_deriv(t) * epsilon;
    P_epsilon = expm(-tau_eps * self.cache.H{t});
else
    H_eps = self.cache.H{t} +(epsilon * self.seq.fields_deriv(t, c)) * self.system.B{c};
    P_epsilon = expm(-self.seq.tau(t) * H_eps);
end

X_n = self.cache.L{t+1} * (P_epsilon * self.cache.U{t});
E_at_eps_point = normalized_distance(self.system.X_final, X_n, self.system.norm2);
ret = (E_at_eps_point - self.cache.E) / epsilon;
