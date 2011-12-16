function term_reason = search_BFGS(user_options)
global OC;


% default options
problem.options = optimset(...
    'TolX',         1e-8,...
    'TolFun',       1e-8,...
    'DerivativeCheck', 'off',...
    'FinDiffType',  'central',...    
    'GradObj',      'on',...
    'LargeScale',   'off', ...
    'OutputFcn',    @monitor_func,...
    'Display',      'off');

% additional user-defined options
if nargin == 1
    problem.options = optimset(problem.options, user_options);
end

% save a copy
OC.opt.BFGS_options = problem.options;


problem.objective = @goal_and_gradient_function_wrapper;
problem.x0 = control_get(OC.opt.control_mask);
problem.solver = 'fminunc';

% try to minimise objective function to zero
[controlFinal, costFinal, exitflag, output] = fminunc(problem.objective, problem.x0, problem.options);

control_update(controlFinal, OC.opt.control_mask); % It may be different than the last point evaluated, and there is no problem-space 

term_reason = OC.opt.term_reason;
end


function [v, grad] = goal_and_gradient_function_wrapper(x)
% Note: x is a vector containing a (possible) subset of the controls
    
    global OC;
    
    OC.opt.N_eval = OC.opt.N_eval + 1;

    control_update(x, OC.opt.control_mask);
    [v, grad] = OC.config.error_func(OC.opt.control_mask);
    OC.opt.last_grad_norm = sum(sum(grad.*grad));
end

