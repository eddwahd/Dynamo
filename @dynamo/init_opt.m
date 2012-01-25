function init_opt(self, control_mask)
% Initialize the optimization module.
% FIXME interplay with the search_* functions and their options, monitor_func.

self.opt.initial_controls = self.seq.get();
self.opt.control_mask = control_mask;


fprintf('Optimization space dimension: %d\n', sum(sum(control_mask)));

self.opt.N_iter = 0;
self.opt.N_eval = 0;
self.opt.last_grad_norm = NaN;
self.opt.term_reason = 'none yet';


self.opt.term_cond = struct( ...
    'max_loop_count',     1e10, ...
    'error_goal',        1e-4, ...
    'max_wall_time',       300, ...
    'max_cputime',         1e4, ...
    'min_gradient_norm', 1e-20);


self.opt.wall_start = now();
self.opt.cpu_start = cputime();


self.stats.error = [];
self.stats.wall_time = [];
self.stats.cpu_time  = [];
self.stats.integral  = [];
self.stats.fluence   = [];
end
