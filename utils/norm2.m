function ret = norm2(A)
% Squared Frobenius norm of an operator (or a vector).

% Ville Bergholm 2011


% real() is just taking care of rounding errors
ret = real(inprod(A, A));
