function [err, grad] = error_open(self, control_mask)
% Error function and its gradient for open (Markovian) systems.
%
% If no control_mask is given, computes just the error function.
% Otherwise also gives the corresponding gradient.

% D(A,B) = d^2(A,B)/|A|^2 = 1 +|B|^2/|A|^2 -2 f(A,B) 


err = normalized_distance(self.system.X_final, self.X(), self.system.norm2);
self.cache.E = err;

if nargin == 2
    grad = self.gradient(control_mask);
end

end
