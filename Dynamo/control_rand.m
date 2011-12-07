function mask = control_rand(T, timeslots, optimize_tau)
% Generates a reasonable set of random initial controls.

global OC;

% shape vectors
normal_controls = [timeslots, length(OC.system.B)];
t_controls = [timeslots, 1];

% Generate random initial controls
initial_controls = 1 * randn(normal_controls);

if optimize_tau
    % optimize tau
    tau_par = [0.25 * T, 0.75 * T]; % min, delta
    %tau_c = 2*pi*rand(t_controls); % random
    tau_c = pi/2 * ones(t_controls); % halfway
    mask = [true(normal_controls), true(t_controls)];
else
    % fixed tau
    tau_par = [T, 0];
    tau_c = zeros(t_controls);
    mask = [true(normal_controls), false(t_controls)];
end

dynamo_init_controls([initial_controls, tau_c], tau_par);
