function [Q, grad_Q] = Q_real(control_mask)
% Nonprojective goal function and gradient.
%
% If no control_mask is given, computes just the goal function.
% Otherwise also gives the corresponding gradient.

global OC;

if nargin == 1
    grad_g = OC.config.gradientFunc(control_mask);    

    grad_Q = real(grad_g) / OC.system.norm2;
end

Q = real(g_func()) / OC.system.norm2;
