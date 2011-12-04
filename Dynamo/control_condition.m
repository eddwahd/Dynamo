function x = control_condition(f)
% Finds the fluence corresponding to the least-ill conditioning for a problem instance.

% Ville Bergholm 2011
  
  global OC;
  q = 10;
  f = logspace(0, 3, 40);
  cond = zeros(size(f));
  for k = 1:length(f)
    cond(k) = find_cond(OC.seq, f(k), q);
  end
  figure
  loglog(f, cond, '.');
  xlabel('fluence')
  ylabel('conditioning')

  % we may assume that cond(fluence) is unimodal
  
  options = optimset('MaxIter', 100, 'MaxFunEvals', 70, 'TolX', 1, 'FunValCheck', 'on');
  %'OutputFcn', @monitor_func,...
  %'Display',   'off');
  [x, fval, eflag] = fminbnd(@(x)find_cond(OC.seq, x, 10), 1, 1e3, options)
end


function ret = find_cond(seq, f, q)
% compute the conditioning for random controls with given fluence

% average over several random control sets
  ret = 0;
  for j = 1:q
    raw = rand_controls_with_fluence(seq, f);

    n_timeslots = size(raw, 1);
    % only update the control fields, not taus
    control_update(raw, [true(size(raw)), false(n_timeslots, 1)]);

    [Q, J] = Q_nr(true(size(seq.raw_controls)));
    cond = conditioning(Q, J);

    ret = ret + cond;
  end
  ret = ret / q;
end


function a = rand_controls_with_fluence(seq, f)
% Generates a random raw control parameter vector
% (without tau) with fluence f.
% Assumes we may modify all controls.
    
    function r = rand_vector(n)
    % Random vector with |r| = 1.
    % TODO Haar measure
        r = randn(n);
        r = r / norm(r, 'fro');
    end

    s = size(seq.control);
    a = rand_vector(s);

    % M is strictly positive so we can do this
    temp = inv(sqrtm(seq.M));
    
    % multiply every timeslot
    a = sqrt(f) * a * temp;
    for k = 1:s(1)
        a(k,:) = a(k,:) / sqrt(seq.tau(k));
    end
    
    % transformed controls: inverse transform
    n_controls = length(seq.control_type);
    for k = 1:n_controls
      switch seq.control_type(k)
        case '.'  % no transformation
        
        case 'p'  % strictly nonnegative, u_k = r_k^2
          a(:, k) = sqrt(abs(a(:, k)));
        
        case 'm'  % minimum and delta, u_k = min + delta * 0.5 * (1 + cos(r_k))
          error('Not yet implemented.')
      
        otherwise
          error('Unknown control type.')
      end
    end
end




function c = conditioning(L, J)
% Compute the ill-conditioning.
% The norm of the least-squares solution of L + J*p = 0
% (using the left pseudoinverse of J since J has more rows than cols)
%p = (J' * J) \ (J' * -L);
    p = pinv(J);
    c = norm(p);
end


