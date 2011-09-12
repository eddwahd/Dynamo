function v = vec(rho)
% VEC  Flattens a matrix into a vector.
%  v = vec(rho)
%
%  Matrix rho is flattened columnwise into a column vector v.
%
%  Used e.g. to convert state operators to superoperator representation.

% JDW 2009
% Ville Bergholm 2009

v = rho(:);
