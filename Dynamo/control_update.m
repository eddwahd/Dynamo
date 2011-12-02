function control_update(x, control_mask)
% control_update  Updates selected controls.
%
%  x: vector of control values, control_mask: corresponding mask.
%
%  Updates control values for which control_mask is true.
%  Makes the changed timeslots and stuff that depends on them stale.

global OC;

if nargin == 1
    control_mask = true(size(OC.seq.raw_controls)); % full mask
end

% make a trial copy of the new controls
new_controls = OC.seq.raw_controls;
new_controls(control_mask) = x;

% see which timeslots have changed
changed_t_mask = any(new_controls ~= OC.seq.raw_controls, 2);

if any(changed_t_mask)
    % actually update the controls
    OC.seq.raw_controls = new_controls;
    
    OC.cache.H_is_stale(changed_t_mask) = true;
    OC.cache.P_is_stale(changed_t_mask) = true;
    
    % Propagate the H_is_stale to the U and Ls.
    OC.cache.U_is_stale( (find(OC.cache.H_is_stale, 1, 'first')+1):end) = true;
    OC.cache.L_is_stale(1:find(OC.cache.H_is_stale, 1, 'last'))         = true;
    
    OC.cache.g_is_stale = true;
    OC.cache.g = NaN;

    % transform the controls
    [OC.seq.tau, OC.seq.tau_deriv, OC.seq.control, OC.seq.control_deriv] = control_transform(OC.seq.raw_controls);
end

