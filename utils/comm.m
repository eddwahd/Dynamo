function S = comm(H, q)
% COMM  Superoperator equivalent for a commutator.
%  S = comm(H);
%
%  [H, rho] == inv_vec(comm(H) * vec(rho))

% Ville Bergholm 2011


S = lmul(H) - rmul(H);
