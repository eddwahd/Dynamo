function [err, grad] = error_open(self, control_mask)
% Error function and its gradient for open (Markovian) systems.
%
% If no control_mask is given, computes just the error function.
% Otherwise also gives the corresponding gradient.

% D(A,B) = d^2(A,B)/|A|^2 = 1 +|B|^2/|A|^2 -2 f(A,B) 

if nargin == 2
    temp = self.gradient_open_1st_order(control_mask);
    grad = (2 / self.system.norm2) * temp;
end

X_n = self.g_func(false);

% fidelity
f = real(inprod(self.system.X_final, X_n));

% |X_n|^2
B_norm2 = norm2(X_n);

err = 1 + (B_norm2 -2*f) / self.system.norm2;
end
