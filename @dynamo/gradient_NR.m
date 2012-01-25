function [J] = gradient_NR(self, control_mask, Q)
% Exact Jacobian of the Newton-Raphson quality function L
% with respect to the controls specified in control_mask. Uses eigendecomposition.
% Returns a matrix.


% compute the eigenvalue factors for the target function
[v, zeta] = eig_factors(Q, true);

% request calculations
% (taus and other controls require different things)
slot_mask = any(control_mask, 2);
tau_slot_mask = control_mask(:, end);
c_slot_mask   = any(control_mask(:, 1:end-1), 2);

self.cache.H_needed_now(slot_mask) = true;   % H_{slot}
self.cache.P_needed_now(c_slot_mask) = true; % P_{c_slot}, also gives H_v and H_eig_factor
temp = [c_slot_mask; false] | [false; tau_slot_mask]; % U_{c_slot},  U_{tau_slot+1}
self.cache.U_needed_now(temp) = true;
self.cache.L_needed_now([false; slot_mask]) = true;  % L_{slot+1}

% and perform them
self.cache_refresh();

% tau as last column of controls
tau_c = size(control_mask, 2);

% iterate over the indices of the true elements of control_mask
J = zeros(prod(size(self.system.X_final)), nnz(control_mask));

[Ts, Cs] = ind2sub(size(control_mask), find(control_mask));
for z = 1:length(Ts)
    t = Ts(z);
    c = Cs(z);

    if c == tau_c
        % dP_t/dtau_{t} = -H_t P_t = -P_t H_t
        temp = -self.seq.tau_deriv(t) * self.cache.L{t+1} * self.cache.H{t} * self.cache.U{t+1};
    else
        % Compute the derivative dP_t/du_{tc} using the eigendecomposition of H(t)
        % For efficiency -tau(t) has been factored out and applied later.
        temp = self.cache.H_v{t}' * self.system.B{c} * self.cache.H_v{t};
        dPdu_eigenbasis = temp .* self.cache.H_eig_factor{t}; % note the elementwise multiplication .* here
        dPdu = self.cache.H_v{t} * dPdu_eigenbasis * self.cache.H_v{t}'; % to computational basis

        temp = -self.seq.tau(t) * self.seq.fields_deriv(t, c) * self.cache.L{t+1} * dPdu * self.cache.U{t};
    end
    J(:, z) = vec(project_to_su(v * ((v' * temp * v) ./ zeta) * v'));
end
end
