function [diff, s] = test_gradient(d, seed)
% Checks if the computed gradient of the error function is accurate.
%
% d is a Dynamo instance containing the optimization problem used.
% If no d is given, uses one of the test suite problems.

% Ville Bergholm 2011-2012


if nargin >= 2
    randseed(seed);
end

%% set up a system, random controls

if nargin < 1
    d = test_suite(21);

    %% choose an error function and a compatible gradient

    %d.config.error_func = @error_abs;
    %d.config.error_func = @error_real;
    d.config.error_func = @error_full;
    d.config.dimS = length(d.system.X_initial);

    d.config.epsilon = 1e-3;

    %d.config.gradient_func = @gradient_g_exact;
    %d.config.gradient_func = @gradient_g_1st_order;
    %d.config.gradient_func = @gradient_g_finite_diff;

    %d.config.gradient_func = @gradient_full_1st_order;
    d.config.gradient_func = @gradient_full_finite_diff;
    % TODO explain the asymptotic O(s) behavior of
    % gradient_full_finite_diff at very small delta (epsilon seems
    % to set the cutoff point)
end
mask = d.full_mask(false);


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
legend('gradient error', 'quadratic', 'linear')
end
