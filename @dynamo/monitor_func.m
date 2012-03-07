function stop = monitor_func(self, x, optimValues, state)
% Executed once every iteration during optimization, decides if we should stop here.


stop = false;
wt = (now() - self.opt.wall_start) * 24*60*60; % elapsed time in seconds
ct = cputime() - self.opt.cpu_start;


%% Plot the sequence every now and then

drawnow(); % flush drawing events and callbacks (incl. user interrupts)

% check for user interrupt
if isempty(self.opt.UI_fig)
    self.opt.term_reason = 'User interrupt';
    stop = true;
else
    if self.opt.plot_interval && mod(self.opt.N_iter, self.opt.plot_interval) == 0
        % Plot the current sequence in the UI figure.
        h = self.opt.UI_fig;
        % It's incredible how much work it takes just to make
        % MATLAB not steal window focus when it plots something.
        set(0, 'CurrentFigure', h);
        ax = get(h, 'CurrentAxes');
        set(ax, 'FontSize', 16);
        self.seq.plot(ax, true);
        text(0.05, 0.9, sprintf('Error: %6.6g', optimValues.fval), 'Units','normalized',...
             'FontSize',18, 'BackgroundColor',[.8 .8 1])
        drawnow();
    end
end


%% Check termination conditions

% TODO some of these are already present in optimValues...
self.opt.N_iter = self.opt.N_iter + 1;

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
