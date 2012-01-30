function mask = init_control(self, T, n_timeslots, optimize_tau, const_value, varargin)
% Generates a reasonable set of initial controls.
% If const_value == [], generates random controls.
    
n_controls = length(self.system.B);

% shape vectors
f_shape = [n_timeslots, n_controls];
t_shape = [n_timeslots, 1];

% tau controls
tau_par = [0.5 * T, 1.0 * T]; % min, delta
tau_c = acos(0) * ones(t_shape); % halfway
fprintf('Tau values ');
if optimize_tau
    fprintf('optimized.\n');
    mask = [true(f_shape), true(t_shape)];
else
    fprintf('fixed.\n')
    mask = [true(f_shape), false(t_shape)];
end


%% Create the control sequence

self.seq = control(n_timeslots, n_controls, tau_par, varargin{:});
self.cache_init(n_timeslots);


%% Set initial control values

if ~isempty(const_value)
    % constant initial controls
    if isscalar(const_value)
        initial_controls = const_value * ones(f_shape);
    else
        % row vector, one value for each control field
        initial_controls = ones(n_timeslots, 1) * const_value;
    end
else
    % Generate random initial controls
    initial_controls = randn(f_shape);
end


self.seq = self.seq.set([initial_controls, tau_c]);
end

