function ret = cache_refresh ()
% This function does most of the heavy computing.
% It should be called _after_ setting _all_ the required OC.cache.X_needed_now fields
% to minimize unnecessary computational work.

global OC;

% The 'needed_now' function looks at everything we need to recompute now,
% compared to what in principle is_stale, and (in principle) executes the
% optimal set of operations so that for everything which was marked
% 'needed_now' is up-to-date.
%
% It then updates the 'is_stale' to indicate what has been actually
% computed, and clears the 'needed_now'
%
% Assumption: the inter-dependence of H and U/L updates is taken care of in 'control_update'


U_recompute_now = OC.cache.U_needed_now & OC.cache.U_is_stale;
L_recompute_now = OC.cache.L_needed_now & OC.cache.L_is_stale;

% To recompute U, you need to start at a cell that is fully recomputed.
% U{1} and L{n_timeslots+1} are by definition never stale and never recomputed.
n_timeslots = length(OC.seq.tau);
for t=n_timeslots:(-1):2
    if U_recompute_now(t+1) && OC.cache.U_is_stale(t)
        U_recompute_now(t) = true;
    end
end
for t=2:n_timeslots
    if L_recompute_now(t-1) && OC.cache.L_is_stale(t)
        L_recompute_now(t) = true;
    end
end

% Now that we know which Us & Ls we need to recompute, we can figure out which Ps and Hs must be up-to-date
P_recompute_now = (U_recompute_now(2:end) | L_recompute_now(1:(end-1)) | OC.cache.P_needed_now) & OC.cache.P_is_stale;
H_recompute_now = (P_recompute_now | OC.cache.H_needed_now) & OC.cache.H_is_stale;

% Compute the Hamiltonians
h_idx = find(H_recompute_now);
for t=h_idx
    H = OC.system.A;
    for c = 1:length(OC.system.B)
        u = OC.seq.control(t, c);
        H = H + u * OC.system.B{c};
    end
    OC.cache.H{t} = H;
end

% Compute the exp(H) and any other per-H computation which may be needed for the gradient function
p_idx = find(P_recompute_now);
for t=p_idx
    OC.config.calcPfromHfunc(t); % Compute the Ps - a single piece of propagator
    % Note: calcPfromHfunc may also compute other values which will be needed for gradient calculations
    %       These should be stored in OC.cache. Their up-to-date-ness is identical to that of P.
end

% Compute the Us - forward propagation (we never recompute U{1})
u_idx = find(U_recompute_now);
for t=u_idx
    OC.cache.U{t} = OC.cache.P{t-1} * OC.cache.U{t-1};
end

% Compute the Ls - adjoint system propagation
el_idx = fliplr(find (L_recompute_now));
for t=el_idx
    OC.cache.L{t} = OC.cache.L{t+1} * OC.cache.P{t};
end

% Mark what has been actually computed
temp = [1, length(OC.seq.tau)];

OC.cache.H_is_stale(H_recompute_now) = false;
OC.cache.P_is_stale(P_recompute_now) = false;
OC.cache.U_is_stale(U_recompute_now) = false;
OC.cache.L_is_stale(L_recompute_now) = false;
OC.cache.H_needed_now = false(temp);
OC.cache.P_needed_now = false(temp);
OC.cache.U_needed_now = [false, false(temp)];
OC.cache.L_needed_now = [false(temp), false];

%ret = [sum(H_recompute_now), sum(P_recompute_now), sum(U_recompute_now), sum(L_recompute_now)]

