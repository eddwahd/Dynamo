function [err, grad] = error_open(control_mask)
% Error function and its gradient for open (Markovian) systems.
%
% If no control_mask is given, computes just the error function.
% Otherwise also gives the corresponding gradient.

global OC;


% D(A,B) = d^2(A,B)/|A|^2 = 1 +|B|^2/|A|^2 -2 f(A,B) 

if nargin == 1
    temp = gradient_open_1st_order(control_mask);
    grad = (2 / OC.system.norm2) * temp;
end

% X_n can be computed using any slice k \in [1, n+1]: X_n = L_k * U_k.
% Try to figure out which k requires least additional computation.
k = g_setup_recalc();
cache_refresh();
X_n = OC.cache.L{k} * OC.cache.U{k};

% fidelity
f = real(inprod(OC.system.X_final, X_n));

% |X_n|^2
B_norm2 = norm2(X_n);

err = 1 + (B_norm2 -2*f) / OC.system.norm2;
