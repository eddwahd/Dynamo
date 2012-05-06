function dPdu = dPdu_exact(H_v, H_eig_factor, B)
% Compute the derivative dP_t/du_{tc} using the eigendecomposition of H(t).
% For efficiency -tau(t) has been factored out and applied later.

    temp = H_v' * B * H_v;
    dPdu_eigenbasis = temp .* H_eig_factor; % note the elementwise multiplication .* here
    dPdu = H_v * dPdu_eigenbasis * H_v';    % to computational basis
end
