function [grad] = gradient_exact(control_mask)
% Exact gradient of auxiliary function g.
%
% Gradient of the auxiliary function g with respect to the controls
% specified in control_mask. Uses eigendecomposition.
% Returns a column vector.

global OC;

if ~isfield(OC.cache,'H_v')
    % Prepare for storing the eigendecomposition data
    temp = [1, length(OC.seq.tau)];
    OC.cache.H_v          = cell(temp);
    OC.cache.H_eig_factor = cell(temp);
end

% request calculations
% (taus and other controls require different things)
slot_mask = any(control_mask, 2);
tau_slot_mask = control_mask(:, end);
c_slot_mask   = any(control_mask(:, 1:end-1), 2);

OC.cache.H_needed_now(slot_mask) = true;   % H_{slot}
OC.cache.P_needed_now(c_slot_mask) = true; % P_{c_slot}, also gives H_v and H_eig_factor
temp = [c_slot_mask; false] | [false; tau_slot_mask]; % U_{c_slot},  U_{tau_slot+1}
OC.cache.U_needed_now(temp) = true;
OC.cache.L_needed_now([false; slot_mask]) = true;  % L_{slot+1}

% and perform them
cache_refresh();

% tau as last column of controls
tau_c = size(control_mask, 2);

% iterate over the indices of the true elements of control_mask
grad = zeros(nnz(control_mask), 1);
[Ts, Cs] = ind2sub(size(control_mask), find(control_mask));
for z = 1:length(Ts)
    t = Ts(z);
    c = Cs(z);

    if c == tau_c
        % dP_t/dtau_{t} = -H_t P_t = -P_t H_t
        temp = -OC.seq.tau_deriv(t) * trace_matmul(OC.cache.L{t+1}, OC.cache.H{t} * OC.cache.U{t+1});
    else
        % Compute the derivative dP_t/du_{tc} using the eigendecomposition of H(t)
        % For efficiency -tau(t) has been factored out and applied later.
        temp = OC.cache.H_v{t}' * OC.system.B{c} * OC.cache.H_v{t};
        dPdu_eigenbasis = temp .* OC.cache.H_eig_factor{t}; % note the elementwise multiplication .* here
        dPdu = OC.cache.H_v{t} * dPdu_eigenbasis * OC.cache.H_v{t}'; % to computational basis

        temp = -OC.seq.tau(t) * OC.seq.control_deriv(t, c) * trace_matmul(OC.cache.L{t+1}, dPdu * OC.cache.U{t});
    end
    grad(z) = temp;
end

