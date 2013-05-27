function err = error_real(self, k)
% Nonprojective error function.
% k is the ensemble index.


err = self.config.f_max -real(self.cache.g{k});
self.cache.VUdagger = -1;
