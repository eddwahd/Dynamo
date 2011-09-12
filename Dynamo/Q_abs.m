function [Q, grad_Q] = Q_abs(control_mask)
% Projective goal function and gradient.
%
% If no control_mask is given, computes just the goal function.
% Otherwise also gives the corresponding gradient.

global OC;

if nargin == 1
    % It isn't any cheaper to compute g and grad_g at once rather than sequentially
    % but since g can be computed using any U and L, it might be
    % cheaper to compute the gradient first...
    grad_g = OC.config.gradientFunc(control_mask);
end

g = g_func();
%Q = abs(g)^2 / OC.system.norm2^2;
Q = abs(g) / OC.system.norm2;

if nargin == 1
    % gradient of |g|^2
    %grad_Q = 2 * real(conj(g) * grad_g) / OC.system.norm2^2;

    % gradient of |g|
    grad_Q = real(conj(g)/abs(g) * grad_g) / OC.system.norm2;
end

