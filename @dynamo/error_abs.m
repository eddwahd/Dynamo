function [err, grad] = error_abs(self, control_mask)
% Projective error function and its gradient.
%
% If no control_mask is given, computes just the error function.
% Otherwise also gives the corresponding gradient.

if nargin == 2
    % It isn't any cheaper to compute g and grad_g at once rather than sequentially
    % but since g can be computed using any U and L, it might be
    % cheaper to compute the gradient first...
    grad_g = self.config.gradientFunc(self, control_mask);
end

g = self.g_func();
%err = self.config.f_max^2 - abs(g)^2 / self.system.norm2^2;
err = self.config.f_max - abs(g) / self.system.norm2;

if nargin == 1
    % gradient of |g|^2
    %grad = 2 * real(conj(g) * grad_g) / -self.system.norm2^2;

    % gradient of |g|
    grad = real(conj(g)/abs(g) * grad_g) / -self.system.norm2;
end

end
