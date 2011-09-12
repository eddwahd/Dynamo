function rho = inv_vec(v, dim)
% INV_VEC  Reshapes a vector into a matrix.
%  rho = inv_vec(v) 
%  rho = inv_vec(v, [n m]) 
%
%  Reshapes vector v (length n*m) into a matrix rho (size [n m]),
%  using column-major ordering. If n and m are not given, rho is assumed
%  to be square.
%
%  Used e.g. to convert state operators from superoperator representation
%  to standard matrix representation.

% JDW 2009
% Ville Bergholm 2009


d = length(v);

if (nargin < 2)
  n = sqrt(d);
  if (floor(n) ~= n)
    error('Length of vector v is not a squared integer.')
  end
  dim = [n n];
else
  if (prod(dim) ~= d)
    error('Dimensions n, m are not compatible with given vector v.')
  end
end

rho = reshape(v, dim);
