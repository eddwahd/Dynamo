function err = error_tr(self, k)
% Error function for closed S+E.
% k is the ensemble index.


% trace norm
[U, S, V] = svd(self.cache.g{k});
err = self.config.f_max -trace(S);

% NOTE assumes S is not singular  
self.cache.VUdagger = V * U';
