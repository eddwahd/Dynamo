function ret = gradient_g_exact(self, t, c)
% Exact gradient of the auxiliary function g.

% Gradient of the auxiliary function g with respect to the controls
% specified in control_mask. Uses the eigendecomposition.
%
% Uses H{t}, U{t_c}, U{t_tau+1} and L{t+1}.

if c < 0
    % dP_t/dtau_{t} = -H_t P_t = -P_t H_t
    ret = -self.seq.tau_deriv(t) * trace_matmul(self.cache.L{t+1}, self.cache.H{t} * self.cache.U{t+1});
else
    dPdu = dPdu_exact(self.cache.H_v{t}, self.cache.H_eig_factor{t}, self.system.B{c});
    ret = -self.seq.tau(t) * self.seq.fields_deriv(t, c) * trace_matmul(self.cache.L{t+1}, dPdu * self.cache.U{t});
end
