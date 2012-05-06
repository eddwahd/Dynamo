function ret = gradient_open_1st_order(self, t, c)
% Gradient of error_open by first order approximation.

% dP_k/du_c \approx (-B_c * dt_k) * P_k
% Exact if G_k commutes with B_c.
%
% Uses H{t}, U{t+1} and L{t+1}.

X_n = self.cache.L{t+1} * self.cache.U{t+1};
temp = X_n - self.system.X_final;
    
if c < 0
    % tau control, exact
    ret = -self.seq.tau_deriv(t) * inprod(temp, self.cache.L{t+1} * self.cache.H{t} * self.cache.U{t+1});
else
    % other control, approximate
    self.gradient_test(t, c);
    ret = -self.seq.tau(t) * self.seq.fields_deriv(t, c) * inprod(temp, self.cache.L{t+1} * self.system.B{c} * self.cache.U{t+1});
end
ret = real(ret);
% normalization
ret = (2 / self.system.norm2) * ret;
