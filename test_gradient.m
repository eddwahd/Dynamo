function [err, s] = test_gradient()
% Checks if the computed gradient of the error function is accurate.
% Uses one of the test suite problems.

% Ville Bergholm 2011-2012


% set up a system, random controls
dyn = test_suite(17);
mask = dyn.full_mask(true);

% initial controls
x = dyn.seq.get(mask);

% error function and its gradient at x
[L, J] = dyn.error_NR(mask);

% random direction in parameter space
temp = randn(size(x));
direction = temp / norm(temp);

err = [];
s = logspace(0, -6, 10);
for k=1:length(s)
    delta = s(k) * direction;

    % linear estimate of error func at x+delta
    predicted = L + J * delta;
    
    % error func at x+delta
    dyn.update_controls(x + delta, mask);
    accurate = dyn.error_NR();
    
    err(k) = norm(predicted - accurate);
end

figure();
loglog(s, err, 'b-o'); % error should be O(s^2)
xlabel('|\Delta c|')
ylabel('|error|')
end
