function shake(self, rel_change)
% Makes a small random perturbation to the control sequence, can be used
% to shake the system out of a local optimum.
% Does not touch the tau values.

n_timeslots = self.seq.n_timeslots();
n_controls = self.seq.n_controls();

% shape vectors
f_shape = [n_timeslots, n_controls];

raw = self.seq.get();
raw(:, 1:n_controls) = raw(:, 1:n_controls) .* (1 +rel_change * randn(f_shape));

self.seq.set(raw);
self.cache.invalidate(); % flush the entire cache
end
