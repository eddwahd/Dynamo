function [v, zeta, exp_d] = eig_factors(G, antihermitian)
% Compute the matrix of eigenvectors v and eigenvalue zeta factors of G.
%  [v, zeta] = eig_factors(G)
%
%  Let G = A +\sum_k u_k B_k,  [v, lambda] = eig(G). Now
%
%    <v(a)| \frac{\partial exp(G)}{\partial u_k} |v(b)> == <v(a)| B_k |v(b)> zeta(a, b),
%
%  where zeta(a, b) := (exp(lambda(a)) - exp(lambda(b))) / (lambda(a) - lambda(b)).
%  When lambda(a) = lambda(b), zeta(a, b) := exp(lambda(a)).
%
% References: 
%     [1] arXiv 1011.4874 (http://arxiv.org/abs/1011.4874)
%     [2] T. Levante, T. Bremi, and R. R. Ernst, J. Magn. Reson. Ser. A 121, 167 (1996)
%     [3] K. Aizu, J. Math. Phys. 4, 762 (1963)
%     [4] R. M. Wilcox, J. Math. Phys. 8, 962 (1967)


% NOTE: If a matrix A has degenerate eigenvalues, the MATLAB function eig(A)
% does not usually return orthogonal eigenvectors for the
% degenerate eigenspace, _unless_ A is Hermitian.
    
% HACK: we only handle Hermitian and antihermitian matrices.
% antihermitian ones we multiply by 1i to make them Hermitian
% before using eig...
if antihermitian
    [v, d] = eig(1i * G);
    d = -1i * diag(d);
else
    [v, d] = eig(G);
    d = diag(d);
end

ooo = ones(1, length(d));

temp = d * ooo;
diff = temp - temp.'; % diff(j,k) = lambda(j) - lambda(k)

degenerate_mask = abs(diff) < 1e-10;
diff(degenerate_mask) = 1; % To prevent division by zero in next step

exp_d = exp(d);
temp = exp_d * ooo; % temp(j,k) = exp(lambda(j))
exp_diff = temp - temp.'; % exp_diff(j,k) = exp(lambda(j)) - exp(lambda(k))
zeta = exp_diff ./ diff;
zeta(degenerate_mask) = temp(degenerate_mask); % For degenerate eigenvalues, the factor is just the exponent
