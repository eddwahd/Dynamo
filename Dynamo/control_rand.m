function mask = control_rand(T, timeslots, optimize_tau, const_value)
% Generates a reasonable set of random initial controls.

global OC;

scale = 1; % TODO

% shape vectors
normal_controls = [timeslots, length(OC.system.B)];
t_controls = [timeslots, 1];

if nargin == 4
    % constant initial controls
    if isscalar(const_value)
        initial_controls = scale * const_value * ones(normal_controls);
    else
        % row vector, one value for each control field
        initial_controls = scale * ones(timeslots, 1) * const_value;
    end
else
    % Generate random initial controls
    initial_controls = scale * randn(normal_controls);
end

fprintf('Tau values ');
if optimize_tau
    fprintf('optimized.\n');
    tau_par = [0.25 * T, 0.75 * T]; % min, delta
    tau_c = pi/2 * ones(t_controls); % halfway
    mask = [true(normal_controls), true(t_controls)];
else
    fprintf('fixed.\n')
    tau_par = [T, 0];
    tau_c = zeros(t_controls);
    mask = [true(normal_controls), false(t_controls)];
end

dynamo_init_controls([initial_controls, tau_c], tau_par);
