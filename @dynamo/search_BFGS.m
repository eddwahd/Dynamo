function term_reason = search_BFGS(self, user_options)
% BFGS optimization.

fprintf('\nOptimizing algorithm: BFGS. Running...\n\n'); drawnow;
    
% default options
problem.options = optimset(...
    'MaxIter',      1e4,...
    'TolX',         1e-8,...
    'TolFun',       1e-8,...
    'DerivativeCheck', 'off',...
    'FinDiffType',  'central',...    
    'GradObj',      'on',...
    'LargeScale',   'off', ...
    'OutputFcn',    @(x, optimValues, state) monitor_func(self, x, optimValues, state),...
    'Display',      'off');

% additional user-defined options
if nargin == 2
    problem.options = optimset(problem.options, user_options);
end

% save a copy
self.opt.BFGS_options = problem.options;


problem.objective = @(x) goal_and_gradient_function_wrapper(self, x);
problem.x0 = self.seq.get(self.opt.control_mask);
problem.solver = 'fminunc';

% try to minimise objective function to zero
[x, cost, exitflag, output] = fminunc(problem.objective, problem.x0, problem.options);

self.update_controls(x, self.opt.control_mask); % It may be different than the last point evaluated, and there is no problem-space 

term_reason = self.opt.term_reason;
end


function [v, grad] = goal_and_gradient_function_wrapper(self, x)
% x is a vector containing (a subset of) the controls
    
    self.opt.N_eval = self.opt.N_eval + 1;

    self.update_controls(x, self.opt.control_mask);
    [v, grad] = self.config.error_func(self, self.opt.control_mask);
    self.opt.last_grad_norm = sum(sum(grad .* grad));
end
