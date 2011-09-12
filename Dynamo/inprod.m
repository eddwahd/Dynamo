function ret = inprod(A,B)
% Hilbert-Schmidt inner product of operators (or vectors).

% Ville Bergholm 2011

ret = sum(sum(conj(A) .* B));
% == trace(A'*B), except computed in O(n^2) time as opposed to O(n^3).
