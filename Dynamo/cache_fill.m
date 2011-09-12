function cache_fill()
% Will invalidate everything, then re-calc everything in the cache
% Used mostly for debugging (since it essentially overrides all matrix-op optimization mechanisms)

global OC;

cache_invalidate();
OC.cache.H_needed_now(:) = true;
OC.cache.P_needed_now(:) = true;
OC.cache.U_needed_now(:) = true;
OC.cache.L_needed_now(:) = true;

cache_refresh();
g_func();
