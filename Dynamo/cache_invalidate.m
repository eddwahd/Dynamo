function cache_invalidate()
% Invalidates the entire cache.

global OC;

OC.cache.H_is_stale(:) = true;
OC.cache.P_is_stale(:) = true;
OC.cache.U_is_stale(2:end) = true;
OC.cache.L_is_stale(1:(end-1)) = true;
OC.cache.g_is_stale = true;
