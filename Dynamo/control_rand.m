function mask = control_rand(T, timeslots, optimize_tau, const)
% Generates a reasonable set of random initial controls.

global OC;

if nargin < 4
    const = false;
end

scale = 1; % TODO

% shape vectors
normal_controls = [timeslots, length(OC.system.B)];
t_controls = [timeslots, 1];

if const
    initial_controls = scale * ones(normal_controls);
else
    % Generate random initial controls
    initial_controls = scale * randn(normal_controls);
end

if optimize_tau
    disp('Tau values optimized.')
    tau_par = [0.25 * T, 0.75 * T]; % min, delta
    %tau_c = 2*pi*rand(t_controls); % random
    tau_c = pi/2 * ones(t_controls); % halfway
    mask = [true(normal_controls), true(t_controls)];
else
    disp('Fixed tau values.')
    tau_par = [T, 0];
    tau_c = zeros(t_controls);
    mask = [true(normal_controls), false(t_controls)];
end

dynamo_init_controls([initial_controls, tau_c], tau_par);
