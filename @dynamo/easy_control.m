function easy_control(self, fields, randomness, t_dependent)
% Generates a reasonable set of initial controls.
% The control fields are set to the given constant values, with an
% optional random perturbation.

n_timeslots = self.seq.n_timeslots();
n_controls = self.seq.n_controls();

% shape vectors
f_shape = [n_timeslots, n_controls];
t_shape = [n_timeslots, 1];

if nargin < 4
    t_dependent = false;
end
if nargin < 3
    randomness = 0;
end


%% Set initial control values

% tau controls
tau_c = acos(0) * ones(t_shape); % halfway

if ~isempty(fields)
    % constant initial controls
    % row vector, one value for each control field
    
    raw = self.seq.inv_transform(fields);
    if ~t_dependent
        % perturb each control field randomly
        raw = raw .* (1 + randomness * randn(1, n_controls));
    end
    initial_controls = ones(t_shape) * raw;
    if t_dependent
        % add some time-dependent noise on top of the raw controls
         initial_controls = initial_controls .* (1 + randomness * randn(f_shape));
    end
else
    % Generate totally random initial controls
    initial_controls = randn(f_shape);
end

self.seq = self.seq.set([initial_controls, tau_c]);
end
