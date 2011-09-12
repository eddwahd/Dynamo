function S = lmul(L, q)
% LMUL  Superoperator equivalent for multiplying from the left.
%  S = lmul(L);
%
%  L*rho == inv_vec(lmul(L)*vec(rho))

% Ville Bergholm 2009


if (nargin < 2)
  q = size(L, 2); % assume target is a square matrix
end

S = kron(speye(q), L);
