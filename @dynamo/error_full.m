function [err, grad] = error_full(self, control_mask)
% Error function and its gradient for open (Markovian) systems.
%
% If no control_mask is given, computes just the error function.
% Otherwise also gives the corresponding gradient.


X_S = partial_trace(self.X(), self.system.dimSE, 2);
err = 0.5 * norm2(self.system.X_final -X_S) / self.system.norm2;
self.cache.E = err;

if nargin == 2
    grad = self.gradient(control_mask) / self.system.norm2;
end
end
