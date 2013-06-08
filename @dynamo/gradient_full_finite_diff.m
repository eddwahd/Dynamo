function ret = gradient_full_finite_diff(self, t, k, c)
% Gradient of error_full by finite difference method.

% f'(x) = (f(x + eps) -f(x))/eps
% Trivial and relatively slow, but a good reference point.
%
% Uses H{t}, U{t} and L{t+1}.


[P, epsilon] = self.finite_diff_P(t, k, c);

% We may compute the finite diff approximation at three locations:
if 0
    % at P
    dPdu = (P -self.cache.P{t, k}) / epsilon;
    dX_Sdu = partial_trace(self.cache.L{t+1, k} * (dPdu * self.cache.U{t, k}), self.system.dimSE, 2);
else
    X_S = partial_trace(self.cache.L{t+1, k} * (P * self.cache.U{t, k}), self.system.dimSE, 2);    
    if 0
        % at X_S
        dX_Sdu = (X_S -self.cache.g{k}) / epsilon;
    else
        % at E (this seems to be the best choice)
        E = 0.5 * norm2(X_S -self.system.X_final);
        ret = (E -self.cache.E) / epsilon;
        return
    end
end
ret = inprod((self.cache.g{k} -self.system.X_final), dX_Sdu);
