function [grad] = gradient_open_finite_diff(self, control_mask, epsilon)
% Gradient of error_open by finite difference method.
%
% E_open'(x) = (E_open(x + eps) - E_open(x))/eps
% Trivial and relatively slow, but a good reference point.

% Mask which Hs, Us, & Ls we need for this calculation
% (it's more efficient to do so before we ask for the current_value, since then get_current_value's call to cache_refresh
% will be most efficient as it knows of all calculations needed at once, and not piece-meal).

slot_mask = any(control_mask, 2);
self.cache.H_needed_now(slot_mask) = true;           % H_{slot}
self.cache.U_needed_now([slot_mask; false]) = true;  % U_{slot}
self.cache.L_needed_now([false; slot_mask]) = true;  % L_{slot+1}
self.cache_refresh();

E = self.cache.E;

% tau as last column of controls
tau_c = size(control_mask, 2);

grad = zeros(nnz(control_mask), 1);
[Ts, Cs] = ind2sub(size(control_mask), find(control_mask));
for z = 1:length(Ts)
    t = Ts(z);
    c = Cs(z);

    % compute everything using just U{t}, L{t+1} and H{t}:
    if c == tau_c
        tau_eps = self.seq.tau(t) +self.seq.tau_deriv(t) * epsilon;
        P_epsilon = expm(-tau_eps * self.cache.H{t});
    else
        H_eps = self.cache.H{t} +(epsilon * self.seq.fields_deriv(t, c)) * self.system.B{c};
        P_epsilon = expm(-self.seq.tau(t) * H_eps);
    end

    X_n = self.cache.L{t+1} * (P_epsilon * self.cache.U{t});
    E_at_eps_point = normalized_distance(self.system.X_final, X_n, self.system.norm2);

    temp = (E_at_eps_point - E) / epsilon;
    grad(z) = temp;
end
end
