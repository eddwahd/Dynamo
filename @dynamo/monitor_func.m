function stop = monitor_func(self, x, optimValues, state)
% Executed once every iteration during optimization, decides if we should stop here.
    

stop = false;

wt = (now() - self.opt.wall_start) * 24*60*60; % elapsed time in seconds
ct = cputime() - self.opt.cpu_start;


% TODO some of these are already present in optimValues...
self.opt.N_iter = self.opt.N_iter + 1;

% plot the sequence every now and then
if self.opt.plot_interval && mod(self.opt.N_iter, self.opt.plot_interval) == 1
  cla();
  self.seq.plot();
  drawnow();
end

%% Check termination conditions

if self.opt.N_eval >= self.opt.term_cond.max_loop_count
    self.opt.term_reason = 'Loop count limit reached';
    stop = true;
end

if wt >= self.opt.term_cond.max_wall_time
    self.opt.term_reason = 'Wall time limit reached';
    stop = true;
end

if ct >= self.opt.term_cond.max_cputime
    self.opt.term_reason = 'CPU time limit reached';
    stop = true;
end

if self.opt.last_grad_norm <= self.opt.term_cond.min_gradient_norm
    self.opt.term_reason = 'Minimal gradient norm reached';
    stop = true;
end

% have we reached our goal?
if optimValues.fval <= self.opt.term_cond.error_goal
    self.opt.term_reason = 'Goal achieved';
    stop = true;
end


%% Stats collector part

self.stats.error(end+1) = optimValues.fval;
self.stats.wall_time(end+1) = wt;
self.stats.cpu_time(end+1)  = ct;
self.stats.integral(end+1,:) = self.seq.integral();
self.stats.fluence(end+1)    = self.seq.fluence(self.system.M);
end
