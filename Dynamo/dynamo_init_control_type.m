function dynamo_init_control_type(control_type, control_par)
% Initialize control types.

global OC;


% Number of controls.
n_controls = length(OC.system.B);


% Check control types
if nargin < 1
    control_type = char(zeros(1, n_controls) + '.');
elseif length(control_type) ~= n_controls
    error('Length of the control_type vector does not match the number of controls.')
end
OC.seq.control_type = control_type;

if any(control_type(OC.system.B_is_superop) == '.')
    disp('Warning: Liouvillian control ops with possibly negative control values.')
end


% Check control parameters
% For now control_par doesn't have separate entries for each timeslot
if nargin < 2
    control_par = {};
end
OC.seq.control_par = control_par;
