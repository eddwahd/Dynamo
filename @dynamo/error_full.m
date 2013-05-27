function err = error_full(self, k)
% Error function for open (Markovian) systems.
% k is the ensemble index.


% system state at t_n
temp = self.cache.U{end, k};

X_S = partial_trace(temp, self.system.dimSE, 2);
err = 0.5 * norm2(self.system.X_final -X_S);
self.cache.E = err;
