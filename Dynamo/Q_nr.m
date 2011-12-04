function [Q, J] = Q_nr(control_mask)
% The Newton-Raphson (vector valued) goal function and its Jacobian.
% Q is denoted by fraktur L in the docs

% Ville Bergholm 2011


global OC
% TODO instead of g, cache Q? now there's no benefit in storing an
% intermediate value or is there?
%if ~OC.cache.g_is_stale 
%  g = OC.cache.g;
%  return;
%end

% g can be computed using any slice k \in [1, n+1]: g = trace(L_k * U_k).
% Try to figure out which k requires least additional computation.
k = g_setup_recalc();
cache_refresh();
Q = P(logm(OC.cache.L{k} * OC.cache.U{k}));
  
%OC.cache.g_is_stale = false;
%OC.cache.g = g;
    
if nargin == 1
    J = gradient_nr(control_mask, Q);
end

Q = vec(Q);
end


function A = P(A)
% Projects A from u(n) into the traceless subalgebra, su(n).
% Essentially this eliminates the global phase.

  n = length(A);
  temp = trace(A) / n;
  A = A - temp * eye(n);
end

