function [tau, tau_deriv, control, control_deriv] = control_transform(raw_control)
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
% Stretchy bins with min and max duration:  t_raw = 0 <=> min duration
tau = OC.seq.tau_par(:,1) +0.5 * OC.seq.tau_par(:,2) .* (1-cos(t_raw));
tau_deriv = 0.5 * OC.seq.tau_par(:,2) .* sin(t_raw);


temp = size(raw_control) - [0, 1];
control = zeros(temp);
control_deriv = zeros(temp);

n_controls = length(OC.seq.control_type);

for k=1:n_controls
    switch OC.seq.control_type(k)
      case '.'  % no transformation
        control(:, k) = raw_control(:, k);
        control_deriv(:, k) = 1;
        
      case 'p'  % strictly nonnegative, u_k = r_k^2
        control(:, k) = raw_control(:, k).^2;
        control_deriv(:, k) = 2 * raw_control(:, k);
        
      case 'm'  % minimum and delta, u_k = min + delta * 0.5 * (1 - cos(r_k))
        par = OC.seq.control_par{k};
        control(:, k) = par(1) +par(2) * 0.5 * (1 - cos(raw_control(:, k)));
        control_deriv(:, k) = par(2) * 0.5 * sin(raw_control(:, k));
      
      otherwise
        error('Unknown control type.')
    end
end


% TODO u_x = A*cos(phi), u_y = A*sin(phi) (not diagonal, but within the same bin, i.e. block diagonal)
% control_deriv(slot, c_j, raw_k) = d c_{sj} / d raw_{sk}
