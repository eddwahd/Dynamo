function ret = gradient_g_finite_diff(self, t, k, c)
% Gradient of the auxiliary function g by finite difference method.
% Special case of gradient_tr_finite_diff.

% g'(x) = (g(x + eps) - g(x))/eps
% Trivial and relatively slow, but a good reference point.
%
% Uses g, H{t}, U{t} and L{t+1}.

epsilon = self.config.epsilon;
if c < 0
    tau_eps = self.seq.tau(t) +self.seq.tau_deriv(t) * epsilon;
    P_epsilon = expm(-tau_eps * self.cache.H{t, k});        
else
    H_eps = self.cache.H{t, k} +(epsilon * self.seq.fields_deriv(t, c)) * self.system.B{c, k};
    P_epsilon = expm(-self.seq.tau(t) * H_eps);
end

g_at_eps_point = trace_matmul(self.cache.L{t+1, k}, P_epsilon * self.cache.U{t, k});
ret = (g_at_eps_point -self.cache.g{k}) / epsilon;

ret = ret * -self.cache.VUdagger;
