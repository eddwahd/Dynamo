function err = error_full(self, k)
% Error function for open (Markovian) systems.
% k is the ensemble index.


temp = self.cache.g{k} -self.system.X_final;
err = 0.5 * norm2(temp);

self.cache.E = err;
