function gradient_test(self, t, c)
% Checks whether the gradient approximation has broken down
% at time slice t, control c.

% FIXME check, what about c?
% The approximation is exact when this vanishes.
violation = abs(self.seq.tau(t)) * norm(self.cache.H{t});

if violation > self.opt.max_violation
    self.opt.max_violation = violation;
    
    if violation > 1    
        fprintf('warning: gradient approximation not valid at t = %d, c = %d: violation = %f.\n', t, c, violation)
    end
    %if violation > 1e2
    %error('damn it')
    %end
end
