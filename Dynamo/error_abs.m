function [err, grad] = error_abs(control_mask)
% Projective error function and its gradient.
%
% If no control_mask is given, computes just the error function.
% Otherwise also gives the corresponding gradient.

global OC;

if nargin == 1
    % It isn't any cheaper to compute g and grad_g at once rather than sequentially
    % but since g can be computed using any U and L, it might be
    % cheaper to compute the gradient first...
    grad_g = OC.config.gradientFunc(control_mask);
end

g = g_func();
%err = OC.system.max_Q^2 - abs(g)^2 / OC.system.norm2^2;
err = OC.system.max_Q - abs(g) / OC.system.norm2;

if nargin == 1
    % gradient of |g|^2
    %grad = 2 * real(conj(g) * grad_g) / -OC.system.norm2^2;

    % gradient of |g|
    grad = real(conj(g)/abs(g) * grad_g) / -OC.system.norm2;
end

