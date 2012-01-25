function term_reason = search_NR(self)
% Newton-Raphson search

fprintf('\nOptimizing algorithm: Newton-Raphson. Running...\n\n'); drawnow;
    
optimValues = struct('fval', 'nonsense');

% which controls are we allowed to change?
mask = self.opt.control_mask;


%control_condition(mask); % preconditioned initial controls

r = 9;  % initial guess for trust region size? needs research

for k = 1:50
  [L, J] = error_NR(mask);
  err = norm(L)^2; % error function to minimize

  x = control_get(mask);
  optimValues.fval = err;
  if monitor_func(x, optimValues, NaN)
    disp('done and done')
    break
  end

  % least squares solution within the trust region
  p = convex_solve(J, L, r);

  % linearized model prediction of the error func value at the new point
  model = norm(J*p + L)^2;

  % update controls, compute actual value
  control_update(x + p, mask);
  actual = norm(error_NR())^2;
  
  % update trust region size
  r = r_update(r, model, actual, err);

  self.opt.N_eval = self.opt.N_eval + 1;
  %self.opt.last_grad_norm = sum(sum(grad.*grad));
end
end

       
function r = r_update(r, model, actual, err)
% Updates the trust region size.

  % Pierre's mystery code for updating r
  temp = -(model - err) / abs(actual - err);
  if temp > 1
    r = r * 0.5;
  elseif temp > 0.3
    r = r * 0.75;
  elseif temp < 0.2
    r = r * 1.25;
  end
end


function x = convex_solve(J, L, r)
% Solves the convex optimization problem
%    min_x |J*x + L|^2  subject to |x|^2 <= r.
%    If the solution is degenerate, return x with the minimal norm?
%    See eqs (5) and (6).

  s = size(J);

  %% Using the CVX convex optimization package...
  q = cvx_quiet(true); % shut up
  cvx_begin
    variable x(s(2));
    minimize( norm(J*x + L) );
    subject to
      norm(x) <= r;
  cvx_end
  cvx_quiet(q); % restore quietness state

  %% done.
    
  % check to see if optimal point is orthogonal to ker J (seems to be)
  %Z = null(J);
  %max(abs(x'*Z))
end
