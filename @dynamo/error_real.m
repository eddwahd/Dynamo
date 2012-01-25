function [err, grad] = error_real(self, control_mask)
% Nonprojective error function and its gradient.
%
% If no control_mask is given, computes just the error function.
% Otherwise also gives the corresponding gradient.

if nargin == 2
    grad_g = self.config.gradient_func(self, control_mask);    
    grad = real(grad_g) / -self.system.norm2;
end

err = self.config.f_max - real(self.g_func()) / self.system.norm2;
end
