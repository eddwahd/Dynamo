function [err, grad] = error_real(control_mask)
% Nonprojective error function and its gradient.
%
% If no control_mask is given, computes just the error function.
% Otherwise also gives the corresponding gradient.

global OC;

if nargin == 1
    grad_g = OC.config.gradientFunc(control_mask);    
    grad = real(grad_g) / -OC.system.norm2;
end

err = OC.system.f_max - real(g_func()) / OC.system.norm2;
