function ret = control_fluence()
% Computes the total fluence of the control fields corresponding to control parameters a.
% The fluence is defined as F^2 = \int_0^T |H_c(t)|^2 dt,    H_c(t) = \sum_r H_r f_r(t).
% If the control Hamiltonians are orthonormal, (H_i, H_j) = \delta_{ij},
% this gives F^2 = \sum_k \int_0^T |f_k(t)|^2 dt.
%
% For piecewise constant controls this reduces to F^2 = \sum_{ri} |a_{ri}|^2 \Delta t_{i}
% Fluence has units of sqrt(energy^2 * time).

% Ville Bergholm 2011

global OC;



ret = 0;
n_timeslots = size(OC.seq.control, 1);
for k = 1:n_timeslots
    c = OC.seq.controls(k, :); % all controls for this timeslot
    ret = ret +c' * OC.seq.M * c * OC.seq.tau{k};
end
