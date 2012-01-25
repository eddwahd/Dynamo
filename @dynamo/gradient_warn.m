function gradient_warn(self, t, c, violation)
% Prints a warning that the gradient approximation has broken down
% at time slice t, control c.


if violation > self.opt.max_violation
    fprintf('warning: gradient approximation not valid at t = %d, c = %d: violation = %f.\n', t, c, violation)
    self.opt.max_violation = violation;
end
if violation > 1e2
    %error('damn it')
end
