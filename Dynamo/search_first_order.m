function term_reason = search_first_order()
global OC;

if ~isfield(OC.config, 'FirstOrder')
    OC.config.FirstOrder = struct();
end
if ~isfield(OC.config.FirstOrder,'step_size')
    OC.config.FirstOrder.step_size = 0.1; 
end

stop = false;

x = controls_get(OC.opt.control_mask);

while ~stop
    OC.opt.N_eval = OC.opt.N_eval + 1;

    [v, grad] = OC.config.Q_func(OC.opt.control_mask);
    
    x = x + OC.config.FirstOrder.step_size .* grad(:);
    OC.opt.last_grad_norm = sqrt(sum(sum(grad .* grad)));
    controls_update(x, OC.opt.control_mask);

    next_v = OC.config.Q_func();
    
    actual_improvement = next_v - v;
    exptected_improvement = OC.opt.last_grad_norm.^2 * OC.config.FirstOrder.step_size / OC.system.norm2;
    
    if actual_improvement < (4/12) * exptected_improvement
        OC.config.FirstOrder.step_size = OC.config.FirstOrder.step_size * 0.99;
    elseif actual_improvement > (8/12) * exptected_improvement
        OC.config.FirstOrder.step_size = OC.config.FirstOrder.step_size * 1.01;
    end

    stop = monitor_func(x);
end
