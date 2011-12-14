function [grad] = gradient_open_1st_order(control_mask)
% Gradient of error_open by first order approximation.

% dP_k/d_u \approx (-B_u * dt_k) * P_k
% Exact if G_k commutes with B_u.

global OC;

% Mask which Hs, Us, & Ls we need for this calculation
% (it's more efficient to do so before we ask for the current_value, since then get_current_value's call to cache_refresh
% will be most efficient as it knows of all calculations needed at once, and not piece-meal).

slot_mask = any(control_mask, 2);

OC.cache.H_needed_now(slot_mask) = true; % H_{slot}
OC.cache.U_needed_now([false; slot_mask]) = true;  % U_{slot+1}
OC.cache.L_needed_now([false; slot_mask]) = true;  % L_{slot+1}

cache_refresh();

% tau as last column of controls
tau_c = size(control_mask, 2);

grad = zeros(nnz(control_mask), 1);
[Ts, Cs] = ind2sub(size(control_mask), find(control_mask));
for z = 1:length(Ts)
    t = Ts(z);
    c = Cs(z);

    % compute everything using just U{t+1}, L{t+1} and H{t}:

    temp = abs(OC.seq.tau(t)) * norm(OC.cache.H{t});
    if temp > 1
        if temp > OC.opt.max_violation
            fprintf('warning: gradient approximation not valid at t = %d, c = %d: temp = %f.\n', t, c, temp)
            OC.opt.max_violation = temp;
        end
        if temp > 1e2
            error('damn it')
        end
    end

    X_n = OC.cache.L{t+1} * OC.cache.U{t+1};
    temp = X_n - OC.system.X_final;
    
    if c == tau_c
        temp = -OC.seq.tau_deriv(t) * inprod(temp, OC.cache.L{t+1} * OC.cache.H{t} * OC.cache.U{t+1});
    else
        temp = -OC.seq.tau(t) * OC.seq.control_deriv(t, c) * inprod(temp, OC.cache.L{t+1} * OC.system.B{c} * OC.cache.U{t+1});
    end
    grad(z) = real(temp);
end

