function ret = normalized_distance(A, B, norm2A)
% Returns the distance measure based on the squared, normalized Frobenius norm.
%
%   D(A, B) = |A - B|_F^2 / |A|_F^2.

% Ville Bergholm 2012

% Normally norm2(A) is precomputed so we can save a few cpu cycles here.
ret = 1 +(norm2(B) -2 * real(inprod(A, B))) / norm2A;
