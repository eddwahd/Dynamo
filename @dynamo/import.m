function import(self, A)
% Imports a control sequence array into a Dynamo instance.
%    
% The rows of A correspond to pulses, each row consists of
% [tau (s), f_rabi_1 (Hz), phi_1, f_rabi_2 (Hz), phi_2, ...]

% Ville Bergholm 2013

% our internal time unit (in s)
TU = self.system.TU;

% split it up
tau = A(:, 1) / TU;
T = sum(tau)

A = A(:, 2:end);
n_bins = size(A, 1)
n_controls = size(A, 2)/2


f = zeros(size(A));
% transform the controls
for k=1:2:2*n_controls
    amp = A(:, k) * TU * 2*pi; 
    phi = A(:, k+1);
    f(:, k)   = amp .* cos(phi);
    f(:, k+1) = amp .* sin(phi);
end

self.update_controls([f, tau], [], true);
end
