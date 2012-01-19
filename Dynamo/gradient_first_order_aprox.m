function [grad] = gradient_first_order_aprox (control_mask)
% Gradient of auxiliary function g by first order approximation.

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
    
    temp = abs(OC.seq.tau(t)) * norm(OC.cache.H{t});
    if temp > 1
        gradient_warn(t, c, temp);
    end

    if c == tau_c
        temp = -OC.seq.tau_deriv(t) * trace_matmul(OC.cache.L{t+1}, OC.cache.H{t} * OC.cache.U{t+1});
    else
        temp = -OC.seq.tau(t) * OC.seq.control_deriv(t, c) * trace_matmul(OC.cache.L{t+1}, OC.system.B{c} * OC.cache.U{t+1});
    end
    grad(z) = temp;
end

