function [out, desc] = export(self)
% Export a control sequence from Dynamo into an array.
%
% Converts the control fields to polar coordinates (Rabi amplitude, phase).
% The output array rows correspond to pulses/bins, each row consists of 
% [tau (s), f_rabi_1 (Hz), \phi_1, f_rabi_2 (Hz), \phi_2, ...]

% Ville Bergholm 2013

desc = sprintf(['%s\n'...
                'Generated: %s UTC\n'...
                'Rows correspond to pulses, each row consists of\n'...
                '[tau (s), f_rabi_1 (Hz), phi_1, f_rabi_2 (Hz), phi_2, ...]'],...
                self.system.description, self.config.date_UTC);

f = self.seq.fields;
TU = self.system.TU;

n_bins = size(f, 1);
n_controls = size(f, 2)/2;
for k=1:2:2*n_controls
    x = f(:, k);
    y = f(:, k+1);
    a(:, k) = sqrt(x.^2 +y.^2) / TU / (2*pi);
    a(:, k+1) = atan2(y, x);
end

out = [self.seq.tau * TU, a];
end
