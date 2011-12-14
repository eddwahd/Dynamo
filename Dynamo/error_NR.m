function [L, J] = error_NR(control_mask)
% The Newton-Raphson (vector valued) error function L and its Jacobian.
% L is denoted by fraktur L in the docs.

% Ville Bergholm 2011


global OC
% TODO instead of g, cache L? now there's no benefit in storing an
% intermediate value or is there?
%if ~OC.cache.g_is_stale 
%  g = OC.cache.g;
%  return;
%end

% g can be computed using any slice k \in [1, n+1]: g = trace(L_k * U_k).
% Try to figure out which k requires least additional computation.
k = g_setup_recalc();
cache_refresh();
L = logm(OC.cache.L{k} * OC.cache.U{k});
  
if nargin == 1
    J = gradient_nr(control_mask, L);
end

L = vec(P(L));
%OC.cache.g_is_stale = false;
%OC.cache.g = g;
end


function A = P(A)
% Projects A from u(n) into the traceless subalgebra, su(n).
% Essentially this eliminates the global phase.

  n = length(A);
  temp = trace(A) / n;
  A = A - temp * eye(n);
end
