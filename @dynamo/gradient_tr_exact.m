function ret = gradient_tr_exact(self, t, c)
% Exact gradient of error_tr.

% Uses the eigendecomposition.
%
% Uses H{t}, U{t_c}, U{t_tau+1} and L{t+1}.

if c < 0
    % dP_t/dtau_{t} = -H_t P_t = -P_t H_t
    dQdu = partial_trace(self.cache.L{t+1} * self.cache.H{t} * self.cache.U{t+1}, self.system.dimSE, 1);
    ret = -self.seq.tau_deriv(t) * trace_matmul(self.cache.VUh, dQdu);
else
    dPdu = dPdu_exact(self.cache.H_v{t}, self.cache.H_eig_factor{t}, self.system.B{c});
    dQdu = partial_trace(self.cache.L{t+1} * dPdu * self.cache.U{t}, self.system.dimSE, 1);
    ret = -self.seq.tau(t) * self.seq.fields_deriv(t, c) * trace_matmul(self.cache.VUh, dQdu);
end
