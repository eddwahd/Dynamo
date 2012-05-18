function ret = gradient_full_1st_order(self, t, c)
% Gradient of error_full by first order approximation.

% dP_k/du_c \approx (-B_c * dt_k) * P_k
% Exact if G_k commutes with B_c.
%
% Uses H{t}, U{t+1} and L{t+1}.

X_S = partial_trace(self.cache.L{t+1} * self.cache.U{t+1}, self.config.dimS, 2);
temp = X_S - self.system.X_final;
    
if c < 0
    % tau control, exact
    ret = -self.seq.tau_deriv(t) * inprod(temp, partial_trace(self.cache.L{t+1} * self.cache.H{t} * self.cache.U{t+1}, self.config.dimS, 2));
else
    % other control, approximate
    self.gradient_test(t, c);
    ret = -self.seq.tau(t) * self.seq.fields_deriv(t, c) * inprod(temp, partial_trace(self.cache.L{t+1} * self.system.B{c} * self.cache.U{t+1}, self.config.dimS, 2));
end
ret = real(ret);
