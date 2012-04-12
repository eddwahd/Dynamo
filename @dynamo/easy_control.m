function easy_control(self, fields, rel_rand, abs_rand, t_dependent)
% Generates a reasonable set of initial controls.
% The control fields are set to the given constant values, with an
% optional random perturbation.

n_timeslots = self.seq.n_timeslots();
n_controls = self.seq.n_controls();

% shape vectors
f_shape = [n_timeslots, n_controls];
t_shape = [n_timeslots, 1];
c_shape = [1, n_controls];

if nargin < 5
    t_dependent = false;
end
if nargin < 4
    abs_rand = 0;
end
if nargin < 3
    rel_rand = 0;
end


%% Set initial control values

% tau controls
tau_c = acos(0) * ones(t_shape); % halfway

if ~isempty(fields)
    % constant initial controls

    if isscalar(fields)
        fields = fields * ones(1, n_controls);
    end
    % now fields should be a row vector, one value for each control field
    
    raw = self.seq.inv_transform(fields);
    if ~t_dependent
        % perturb each control field randomly
        raw = raw +rel_rand * raw .* randn(c_shape) +abs_rand * randn(c_shape);
        raw = ones(t_shape) * raw;
    else
        % add some time-dependent noise on top of the raw controls
        raw = ones(t_shape) * raw;
        raw = raw +rel_rand * raw .* randn(f_shape) +abs_rand * randn(f_shape);
    end
else
    % Generate totally random initial controls
    raw = randn(f_shape);
end

self.seq.set([raw, tau_c]);
self.cache.invalidate(); % flush the entire cache
end
