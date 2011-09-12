function g = g_func()
% Computes the auxiliary function g := trace(X_f^\dagger * X(t_n)).
% Used both for the goal function as well as its gradient.

global OC;

if ~OC.cache.g_is_stale 
    g = OC.cache.g;
    return;
end

% g can be computed using any slice k \in [1, n+1]: g = trace(L_k * U_k).
% Try to figure out which k requires least additional computation.
k = g_setup_recalc();
cache_refresh();

g = trace_matmul(OC.cache.L{k}, OC.cache.U{k});

OC.cache.g_is_stale = false;
OC.cache.g = g;
