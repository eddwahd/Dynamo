function ret = gradient_g_finite_diff(self, t, c, epsilon)
% Gradient of the auxiliary function g by finite difference method.

% g'(x) = (g(x + eps) - g(x))/eps
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

g_at_eps_point = trace_matmul(self.cache.L{t+1}, P_epsilon * self.cache.U{t});
ret = (g_at_eps_point - self.cache.g) / epsilon;
