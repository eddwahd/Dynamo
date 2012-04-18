function [diff, s] = test_gradient(seed)
% Checks if the computed gradient of the error function is accurate.
% Uses one of the test suite problems.

% Ville Bergholm 2011-2012

randseed(seed);


%% set up a system, random controls

d = test_suite(21);
mask = d.full_mask(true);


%% choose an error function and a compatible gradient

%d.config.error_func = @error_abs;
%d.config.error_func = @error_real;
d.config.error_func = @error_open;

epsilon = 1e-7;

%d.config.gradient_func = @gradient_g_exact;
%d.config.gradient_func = @gradient_g_1st_order;
%d.config.gradient_func = @(s, m) gradient_g_finite_diff(s, m, epsilon);

%d.config.gradient_func = @gradient_open_1st_order;
d.config.gradient_func = @(s, m) gradient_open_finite_diff(s, m, epsilon);
% TODO explain the O(s) behavior of open_finite_diff


%% test the accuracy

% save the initial controls
x = d.seq.get(mask);

% error function and its gradient at x
[err, grad] = d.config.error_func(d, mask);

% random direction in parameter space
temp = randn(size(x));
direction = temp / norm(temp);

s = logspace(0, -6, 20);
diff = [];
for k=1:length(s)
    delta = s(k) * direction;

    % linear estimate of error func at x+delta
    predicted = err + grad.' * delta;
    
    % error func at x+delta
    d.update_controls(x + delta, mask);
    accurate = d.config.error_func(d);
    
    diff(k) = norm(predicted - accurate);
end

figure();
% gradient error should be \propto s^2 for an exact or
% finite_diff gradient with a sufficiently small epsilon, and
% \propto s for a 1st order gradient approximation.
loglog(s, diff, 'b-o', s, s.^2, 'r', s, s, 'g');
xlabel('|\Delta u|')
ylabel('|gradient error|')
end
