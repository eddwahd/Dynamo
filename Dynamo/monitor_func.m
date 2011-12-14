function stop = monitor_func(x, optimValues, state)
% Executed once every iteration during optimization, decides if we
% should stop here.
    
global OC;

stop = false;

wt = (now() - OC.opt.wall_start) * 24*60*60; % elapsed time in seconds
ct = cputime() - OC.opt.cpu_start;


% TODO some of these are already present in optimValues...
OC.opt.N_iter = OC.opt.N_iter + 1;


if OC.opt.N_eval >= OC.opt.term_cond.max_loop_count
    OC.opt.term_reason = 'Loop count limit reached';
    stop = true;
end

if wt >= OC.opt.term_cond.max_wall_time
    OC.opt.term_reason = 'Wall time limit reached';
    stop = true;
end

if ct >= OC.opt.term_cond.max_cputime
    OC.opt.term_reason = 'CPU time limit reached';
    stop = true;
end

if OC.opt.last_grad_norm <= OC.opt.term_cond.min_gradient_norm
    OC.opt.term_reason = 'Minimal gradient norm reached';
    stop = true;
end

% have we reached our goal?
if optimValues.fval <= OC.opt.term_cond.error_goal
    OC.opt.term_reason = 'Goal achieved';
    stop = true;
end


% stats collector part
OC.stats.error(end+1) = optimValues.fval;
OC.stats.wall_time(end+1) = wt;
OC.stats.cpu_time(end+1)  = ct;
OC.stats.fluence(end+1) = control_fluence(OC.seq, OC.system);
end
