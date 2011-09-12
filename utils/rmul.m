function S = rmul(R, m)
% RMUL  Superoperator equivalent for multiplying from the right.
%  S = rmul(R);
%
%  rho*R == inv_vec(rmul(R)*vec(rho))

% Ville Bergholm 2009


if (nargin < 2)
  m = size(R, 1); % assume target is a square matrix
end

S = kron(R.', speye(m));
