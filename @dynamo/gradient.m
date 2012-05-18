function [grad] = gradient(self, control_mask)
% Auxiliary function for computing the various gradients.


% Find out which Hs, Us, & Ls we need to compute the gradient.
% It's usually more efficient to do all the required calculations at once, and not piece-meal.
slot_mask = any(control_mask, 2);
self.cache.H_needed_now(slot_mask) = true;           % H_{slot}
self.cache.L_needed_now([false; slot_mask]) = true;  % L_{slot+1}

% TODO FIXME kind of a hack, the gradient funcs themselves should keep this information
temp = self.config.gradient_func;
if isequal(temp, @gradient_g_1st_order)...
        || isequal(temp, @gradient_full_1st_order)
    self.cache.U_needed_now([false; slot_mask]) = true;  % U_{slot+1}

elseif isequal(temp, @gradient_g_exact)...
        || isequal(temp, @gradient_g_mixed_exact)...
        || isequal(temp, @gradient_tr_exact)
    % taus and other controls require different things
    tau_slot_mask = control_mask(:, end);
    c_slot_mask   = any(control_mask(:, 1:end-1), 2);
    temp = [c_slot_mask; false] | [false; tau_slot_mask]; 

    % TODO such a minute difference, maybe we could always use the more inclusive expression?
    if isequal(temp, @gradient_g_mixed_exact)
        self.cache.P_needed_now(slot_mask) = true; % P{slot}
    else
        self.cache.P_needed_now(c_slot_mask) = true; % P_{c_slot}, also gives H_v and H_eig_factor
    end
    self.cache.U_needed_now(temp) = true;        % U_{c_slot},  U_{tau_slot+1}
    
elseif isequal(temp, @gradient_g_finite_diff) || isequal(temp, @gradient_full_finite_diff)
    self.cache.U_needed_now([slot_mask; false]) = true;  % U_{slot}

    if isequal(temp, @gradient_g_finite_diff)
        self.cache.g_needed_now = true;
    end

else
    error('Unknown gradient function.')
end
self.cache_refresh();


% tau are the last column of the controls
tau_c = size(control_mask, 2);

% iterate over the true elements of the mask
grad = zeros(nnz(control_mask), 1);
[Ts, Cs] = ind2sub(size(control_mask), find(control_mask));
for z = 1:length(Ts)
    t = Ts(z);
    c = Cs(z);
    if c == tau_c
        c = -1; % denotes a tau control
    end
    
    grad(z) = self.config.gradient_func(self, t, c);
end
