function X = system_get(k)
%  Controlled system at t_k.
%
% Returns X(t_k).

% Ville Bergholm 2011

global OC;

% U{k} is the system at t = sum(tau(1:(k-1))) = t_{k-1}
OC.cache.U_needed_now(k+1) = true;
cache_refresh();

X = OC.cache.U{k+1};
end

