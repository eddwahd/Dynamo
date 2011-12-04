function calcPfromH_exact_gradient(t)
% Computed once per OC.cache.H{t}, i.e. once per time slot, for a specific settings of all controls for that slot
%

global OC;

minus_dt_H = -OC.seq.tau(t) * OC.cache.H{t};

N = length(minus_dt_H);

%% Compute the eigenvalue factors and eigenvectors of -dt*H

[v, zeta, exp_d] = eig_factors(minus_dt_H, true);

OC.cache.H_v{t} = v;
OC.cache.H_eig_factor{t} = zeta;


%% And finally expm(-dt*H) using the eigendecomposition

OC.cache.P{t} = v * diag(exp_d) * v';
