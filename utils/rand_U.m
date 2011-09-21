function U = rand_U(n)
% RAND_U  Random U(n) matrix.
%  U = rand_U(n)
%
%  Returns a random unitary n*n matrix.
%  The matrix is random with respect to the Haar measure.

%! F. Mezzadri, "How to generate random matrices from the classical compact groups", Notices of the AMS 54, 592 (2007). arXiv.org:math-ph/0609050
% Ville Bergholm 2005-2009


% sample the Ginibre ensemble, p(Z(i,j)) == 1/pi * exp(-abs(Z(i,j))^2),
% p(Z) == 1/pi^(n^2) * exp(-trace(Z'*Z))
Z = (randn(n) + i*randn(n))/sqrt(2); 

[Q, R] = qr(Z); % QR factorization

% eliminate multivaluedness in Q
P = diag(R); P = P./abs(P); 
U = Q*diag(P);
