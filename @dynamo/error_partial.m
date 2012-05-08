function [err, grad] = error_partial(self, control_mask)
% Error function and its gradient for closed S+E
%
% If no control_mask is given, computes just the error function.
% Otherwise also gives the corresponding gradient.


  k = self.cache.g_setup_recalc();
  self.cache_refresh();
  temp = self.cache.L{k} * self.cache.U{k};
  Q = partial_trace(temp, self.config.dimS);

  % trace norm
  [U, S, V] = svd(Q);
  err = 1 -trace(S) / self.system.norm2;

  if nargin == 2
      % NOTE assumes S is not singular  
      self.cache.VUh = V * U';
      grad = real(self.gradient(control_mask)) / -self.system.norm2;
  end
end
