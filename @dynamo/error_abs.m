function err = error_abs(self, k)
% Projective error function, special case (dim E = 1) of error_tr().
% Separately implemented for performance reasons.
% k is the ensemble index.


g = self.cache.g{k};
temp = abs(g);
err = self.config.f_max -temp;

if temp == 0
    self.cache.VUdagger = 0; % |g| not differentiable at this point
else
    self.cache.VUdagger = conj(g) / temp; % == V * U' for [U, S, V] = svd(g)
end
self.cache.E = err;
