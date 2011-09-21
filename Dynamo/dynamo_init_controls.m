function dynamo_init_controls(controls, par, squared)

%  tau: vector with the min. durations of the time slots, or a scalar denoting the min. total time.
%  controls: initial control values, size == [n_timeslots, n_controls].
%  squared: optional, boolean vector denoting which controls should be squared (nonnegative).
%
%  The initial and final states need to be set before this function is called.


global OC;


% Number of time slices for the piecewise constant controls.
n_timeslots = size(controls, 1);

% Number of control fields.
% The last column are the tau controls.
n_controls  = size(controls, 2) - 1;
if n_controls ~= length(OC.system.B)
    error('Number of control Hamiltonians does not match the number of controls given.')
end

fprintf('Timeslots: %d\nControls: %d\n', n_timeslots, n_controls);


%% Check the validity of the parameters. This code needs to change when controls_transform changes.

% Check squared (nonnegative) control fields
if nargin < 3
    squared = false(1, n_controls);
elseif length(squared) ~= n_controls
    error('Length of the squared vector does not match the number of controls.')
end
OC.seq.squared = squared;

if ~all(squared(OC.system.B_is_superop))
    disp('Warning: Liouvillian control ops with possibly negative control values.')
end


% Is par a total time or a vector of delta_t:s?
len_params = size(par, 1);
if len_params == 1
    % divide the parameters uniformly into timeslots
    OC.seq.par = ones(n_timeslots,1) * (par / n_timeslots);
else
    if len_params ~= n_timeslots
        error('Length of par does not match the number of control timeslots given.')
    end
    OC.seq.par = par;
end


%% Set up controls

OC.config.initial_controls = controls;
OC.seq.raw_controls = controls;

% transform the controls
[OC.seq.tau, OC.seq.tau_deriv, OC.seq.control, OC.seq.control_deriv] = controls_transform(controls);


%% Set up caching

% TODO technically we could check here if
% (a) cache values already exists and
% (b) some of them could still be valid...

temp = [1, length(OC.seq.tau)];

% Generator for a time slice, H
OC.cache.H = cell(temp);

% Propagator for a single time slice. Roughly P = expm(-dt * H). Computed by OC.config.calcPfromHfunc
OC.cache.P = cell(temp);  


% Forward and backward propagators U, L
% U{k+1} = P{k} * U{k};
% U{k} is the system at t = sum(tau(1:(k-1))) = t_{k-1}
OC.cache.U = cell(temp + [0, 1]);
% L{k-1} = L{k} * P{k-1};
% L{k} is the adjoint system at t = sum(tau(1:(k-1))) = t_{k-1}
OC.cache.L = cell(temp + [0, 1]);

% U: X_initial propagated forward up to a time instant.
OC.cache.U{1} = OC.system.X_initial;
switch OC.config.task
  case {'task5', 'task6'}
    % NOTE: in these tasks, L is a pure propagator (even though storing them requires more memory).
    OC.cache.L{end} = eye(length(OC.system.X_final));
    
  otherwise
    % L: X_final propagated backward 
    OC.cache.L{end} = OC.system.X_final';
end


% Keep track of what needs re-computation if we want a complete update of everything.
% U{1} and L{end} are never stale and never recomputed.
OC.cache.H_is_stale = true(temp);
OC.cache.P_is_stale = true(temp);
OC.cache.U_is_stale = [false, true(temp)]; % Updates for H via 'controls_update' get propagated automatically
OC.cache.L_is_stale = [true(temp), false];

% Here we indicate which values we need to have up-to-date
% The 'needed_now' function looks at everything we need to recompute now, 
% compared to what in principle is_stale, and (in principle) executes the
% optimal set of operations so that for everything which was marked
% 'needed_now' is up-to-date

% Example: We modified H{3}, so in theory we need to recompute U{4:end}.
% But for our immediate needs we only want U{7} and L{7},
% so we mark U{4:end} as "is_stale", but only 4:7 as "needed_now".

OC.cache.H_needed_now = false(temp);   
OC.cache.P_needed_now = false(temp);   
OC.cache.U_needed_now = [false, false(temp)];
OC.cache.L_needed_now = [false(temp), false]; 

% Note - certain gradient methods may cache additional computations

OC.cache.g_is_stale = true;
OC.cache.g = NaN;

