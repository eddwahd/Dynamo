function S = lrmul(L, R)
% LRMUL  Superoperator equivalent for multiplying both from left and right.
%  S = lrmul(L, R);
%
%  L*rho*R == inv_vec(lrmul(L, R)*vec(rho))

% Ville Bergholm 2009-2010


% L and R fix the shape of rho completely
S = kron(R.', L);
