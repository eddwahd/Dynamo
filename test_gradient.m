function test_gradient()
% Checks if the computed gradient of the error function is accurate.
% Uses one of the test suite problems.

% Ville Bergholm 2011


% set up a system, random controls
control_mask = test_suite(17);
x = control_get(control_mask);

% error function and its gradient at x
[L, J] = error_NR(control_mask);

% random direction in parameter space
direction = randn(size(x));

err = [];
s = logspace(0, -5, 10);
for k=1:length(s)
    delta = s(k) * direction;

    % linear estimate of error func at x+delta
    predicted = L + J * delta;
    
    % error func at x+delta
    control_update(x + delta, control_mask);
    accurate = error_NR();
    
    err(k) = norm(predicted - accurate);
end

figure();
loglog(s, err); % error should be O(s^2)
err
