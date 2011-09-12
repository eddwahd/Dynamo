function [Q, grad_Q] = Q_open(control_mask)
% Goal function and gradient for open (Markovian) systems.
%
% If no control_mask is given, computes just the goal function.
% Otherwise also gives the corresponding gradient.

global OC;


% Q(A,B) = 1 -d^2(A,B)/|A|^2 = 2 f(A,B) - |B|^2/|A|^2 \le 1

if nargin == 1
    temp = gradient_open_1st_order(control_mask);
    grad_Q = (2 / OC.system.norm2) * temp;
end

% X_n can be computed using any slice k \in [1, n+1]: X_n = L_k * U_k.
% Try to figure out which k requires least additional computation.
k = g_setup_recalc();
cache_refresh();
X_n = OC.cache.L{k} * OC.cache.U{k};

% fidelity
f = real(inprod(OC.system.X_final, X_n));

% |X_n|^2 (real takes care of rounding errors)
B_norm2 = real(inprod(X_n, X_n));


Q = (2*f -B_norm2) / OC.system.norm2;

