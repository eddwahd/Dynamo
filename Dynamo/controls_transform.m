function [tau, tau_deriv, control, control_deriv] = controls_transform(raw_control)
%  Diagonal transformation function for the controls.
%
% Given raw controls r_k, computes the transformed controls c_k(r_k) and their derivatives.
% Uses the (fixed) control parameters defined in OC.seq.par and elsewhere.

% Ville Bergholm 2011

global OC;


% Tau control is the last one.
t_raw = raw_control(:, end);

% Returned tau values should always be positive, and not too large
% if we're using the 1st order gradient approximation.
% Stretchy bins with min and max duration:  t_raw = 0 <=> max duration
tau = OC.seq.par(:,1) +0.5 * OC.seq.par(:,2) .* (1+cos(t_raw));
tau_deriv = -0.5 * OC.seq.par(:,2) .* sin(t_raw);

% Squared control fields.
% squared(k) == true means the control term is u_k^2 B_k (for strictly positive controls)
temp = size(raw_control) - [0, 1];
control = zeros(temp);
control(:, ~OC.seq.squared) = raw_control(:, ~OC.seq.squared);
control(:,  OC.seq.squared) = raw_control(:,  OC.seq.squared).^2;

control_deriv = zeros(temp);
control_deriv(:, ~OC.seq.squared) = 1;
control_deriv(:,  OC.seq.squared) = 2 * raw_control(:, OC.seq.squared);


% TODO u_x = A*cos(phi), u_y = A*sin(phi) (not diagonal, but within the same bin, i.e. block diagonal)
% control_deriv(slot, c_j, raw_k) = d c_{sj} / d raw_{sk}