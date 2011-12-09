function test_suite(p)
% Implementation of the test cases from the paper
% S. Machnes et al., arXiv:1011.4874
% Ville Bergholm 2011



%% Basic definitions: Pauli matrices etc.

SX = [0 1; 1 0];
SY = [0 -1i; 1i 0];
SZ = [1 0; 0 -1];
I = speye(2,2);

CNOT = [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0];


%% Define the problem

fprintf('Test suite problem %d:\n  ', p)
if 1 <= p && p <= 12
  % 1-12: Ising chain
  switch p
    case {1, 2, 3, 4}
      q = 2;
      final = CNOT;

    case {5, 6}
      q = 3;
      final = qft(q);
    
    case {7, 8, 9}
      q = 4;
      final = qft(q);

    case {10, 11, 12}
      q = 5;
      final = qft(q);
  end    
  fprintf('Ising chain, %d qubits, XY control.\n', q)
  dim = 2 * ones(1, q);
  H = heisenberg(dim, 2*[0 0 1]);
  C = control(dim, 'xy');

else
  switch p
  case {13, 14}
    dim = [2 2 2 2];
    fprintf('Completely ZZ-coupled graph, 4 qubits, XY control.\n')
    H_CS = heisenberg(dim, 2*[0 0 1]) +mkron(0.5 * SZ, speye(4,4), SZ);
    H = H_CS +mkron(0.5 * SZ, I, SZ, I) +mkron(0.5 * I, SZ, I, SZ);
    C = control(dim, 'xy');
    final = expm(-1j * pi/2 * H_CS); % target: C_4 cluster state
    
  case {15, 16}
    % NV centers in diamond
    q = 2;
    E = 2*pi * [-134.825, -4.725, 4.275, 135.275]; % MHz
    mu = [1, 1/3.5, 1/1.4, 1/1.8];
    % FIXME controls?
    H = diag(E) +2*pi * 135 * diag([1, 0, 0, -1]);
    final = CNOT;

  case {17, 18}
    q = 5;
    dim = 2 * ones(1, q);
    fprintf('Ising chain with Stark shift, %d qubits, uniform XY control.\n', q)
    H = heisenberg(dim, 2*[0 0 1]) + op_sum(dim, @(k) -SZ*(2+k));
    C = {op_sum(dim, 0.5*SX), op_sum(dim, 0.5*SY)};
    final = qft(q);
    
  case 19
    q = 5;
    dim = 2 * ones(1, q);
    fprintf('Heisenberg chain with bias, %d qubits, Z control.\n', q)
    H = heisenberg(dim, 2*[1 1 1]) + op_sum(dim, -10*SX);
    C = control(dim, 'z');
    final = qft(q);
    
  case {20, 21}
    q = 3 +p -20;
    dim = 2 * ones(1, q);
    fprintf('Heisenberg chain, %d qubits, XY control at one end.\n', q)
    H = heisenberg(dim, 2*[1 1 1]);
    C = control(dim, 'xy', q - 2);
    final = rand_U(prod(dim));

  case {22, 23}
    if p == 22
        d = 13;
    else
        d = 7;
    end
    fprintf('A single spin-%g, J_z, J_x control.\n', (d-1)/2)
    J = angular_momentum(d);
    H = J{3}^2;
    C = {J{3}, J{1}};
    final = rand_U(d);
    
  otherwise
    error('Unknown problem.');
end
end
fprintf('\n')
initial = eye(size(final));

dynamo_init('S gate', initial, final, H, C)



%% Optimization options

timeslots = [30, 40, 128, 64, 120, 140, 128, 128, 64, 300, 300, 64,...
             128, 128, 40, 64, 1000, 1000, 300, 64, 128, 100, 50];
T = [2, 2, 3, 4, 6, 7, 10, 12, 20, 15, 20, 25,...
     7, 12, 2, 5, 125, 150, 30, 15, 40, 15, 5];


% Set up random initial controls.
control_mask = control_rand(T(p), timeslots(p), false);
dynamo_init_opt(control_mask);


%% Now do the actual search

fprintf('\nOptimizing algorithm: GRAPE (BFGS 2nd order update scheme, updating all time slices concurrently).\n\n    Please wait, this may take a while... \n\n'); drawnow;

OC.config.BFGS = struct('fminopt', struct('Display', 'off'));
termination_reason = search_BFGS();

analyze();
end
