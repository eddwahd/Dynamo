function term_reason = search_BFGS()
global OC;


%OC.opt.term_reason = 'BFGS: Minimal gradient norm reached'; % If BFGS cannot find the goal, this is because the gradients are pointing nowhere useful


problem.options = optimset(...
    'TolX',         1e-8,...
    'TolFun',       1e-8,...
    'DerivativeCheck', 'off',...
    'FinDiffType',  'central',...    
    'GradObj',      'on',...
    'LargeScale',   'off', ...
    'OutputFcn',    @monitor_func,...
    'Display',      'off');

if isfield(OC.config,'BFGS') && isfield(OC.config.BFGS,'fminopt')
    fn = fieldnames(OC.config.BFGS.fminopt);
    for k=1:length(fn)
        problem.options = optimset(problem.options, fn{k}, OC.config.BFGS.fminopt.(fn{k}));
   end
end    

problem.objective = @goal_and_gradient_function_wrapper;
problem.x0 = control_get(OC.opt.control_mask);
problem.solver = 'fminunc';

% try to minimise objective function to -1
[controlFinal, costFinal, exitflag, output] = fminunc(problem.objective, problem.x0, problem.options);

control_update(controlFinal, OC.opt.control_mask); % It may be different than the last point evaluated, and there is no problem-space 

term_reason = OC.opt.term_reason;
end


function [v, grad] = goal_and_gradient_function_wrapper(x)
% Note: x is a vector containing a (possible) subset of the controls
    
    global OC;
    
    OC.opt.N_eval = OC.opt.N_eval + 1;

    control_update(x, OC.opt.control_mask);
    [v, grad] = OC.config.Q_func(OC.opt.control_mask);
    % fminunc minimizes, we wish to maximize
    v = -v;  grad = -grad;
    OC.opt.last_grad_norm = sum(sum(grad.*grad));
end

