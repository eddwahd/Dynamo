function [slots, T] = test_suite(p)
% Implementation of the test cases from the paper
% S. Machnes et al., arXiv:1011.4874
% Ville Bergholm 2011



%% Basic definitions: Pauli matrices etc.

SX = [0 1; 1 0];
SY = [0 -1i; 1i 0];
SZ = [1 0; 0 -1];
I = speye(2,2);

SP = (SX +1i*SY)/2;

CNOT = [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0];


%% Define the problem

fprintf('Test suite problem %d:\n  ', p)
switch p
  case {1, 2, 3, 4}
    % Ising chain
    q = 2;
    fprintf('Ising chain, %d qubits, XY control.\n', q)
    H = ising(q, 1);
    C = control_XY(q);
    final = CNOT;

  case {5, 6}
    q = 3;
    fprintf('Ising chain, %d qubits, XY control.\n', q)
    H = ising(q, 1);
    C = control_XY(q);
    final = qft(q);
    
  case {7, 8, 9}
    q = 4;
    fprintf('Ising chain, %d qubits, XY control.\n', q)
    H = ising(q, 1);
    C = control_XY(q);
    final = qft(q);

  case {10, 11, 12}
    q = 5;
    fprintf('Ising chain, %d qubits, XY control.\n', q)
    H = ising(q, 1);
    C = control_XY(q);
    final = qft(q);

  case {13, 14}
    q = 4;
    fprintf('Cluster state, %d qubits, XY control.\n', q)
    H_CS = ising(q, 1);
    H_CS = H_CS + mkron(0.5 * SZ, speye(4,4), SZ);
    H = H_CS +mkron(0.5 * SZ, I, SZ, I) +mkron(0.5 * I, SZ, I, SZ);
    C = control_XY(q);
    final = expm(-1j * pi/2 * H_CS);
    
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
    fprintf('Ising chain with Stark shift, %d qubits, uniform XY control.\n', q)
    H = ising(q, 1) + sum_op(q, -SZ, 2+(1:n));
    C{1} = sum_op(q, 0.5*SX);
    C{2} = sum_op(q, 0.5*SY);
    final = qft(q);
    
  case 19
    q = 5;
    fprintf('Heisenberg chain with bias, %d qubits, Z control.\n', q)
    H = heisenberg(q, 1) + sum_op(q, -10*SX);
    C = control_Z(q);
    final = qft(q);
    
  case {20, 21}
    if (p == 20)
        q = 3;
        C = control_XY(q, 1);
    else
        q = 4;
        C = control_XY(q, 2);
    end
    fprintf('Heisenberg chain, %d qubits, XY control at one end.\n', q)
    H = heisenberg(q, 1);
    final = rand_U(2^q);

  case {22, 23}
    if (p == 22)
        d = 13
    else
        d = 7;
    end
    fprintf('Spin-%g chain, Jz/Jx control.\n', (d-1)/2)
    J = angular_momentum(d);
    H = J{3}^2;
    C{1} = J{3};
    C{2} = J{1};
    final = rand_U(d);
    
  otherwise
    error('Unknown problem.');
end
fprintf('\n\n')
initial = eye(size(final));

dynamo_init('task1', initial, final, H, C);


%% Optimization options

slots = [30, 40, 128, 64, 120, 140, 128, 128, 64, 300, 300, 64,...
         128, 128, 40, 64, 1000, 1000, 300, 64, 128, 100, 50];
T = [2, 2, 3, 4, 6, 7, 10, 12, 20, 15, 20, 25,...
     7, 12, 2, 5, 125, 150, 30, 15, 40, 15, 5];

slots = slots(p);
T = T(p);




function [H] = ising(n, J)
% Ising spin chain Hamiltonian for n spin-1/2:s

    temp = J/2 * kron(SZ, SZ);
    N = 2^n;
    H = sparse(N, N);
    for k=2:n
        H = H + mkron(speye(2^(k-2)), temp, speye(2^(n-k)));
    end
end


function [H] = heisenberg(n, J)
% Heisenberg spin chain Hamiltonian for n spin-1/2:s

    temp = J/2 * (kron(SX,SX) + kron(SY,SY) + kron(SZ,SZ));
    N = 2^n;
    H = sparse(N, N);
    for k=2:n
        H = H + mkron(speye(2^(k-2)), temp, speye(2^(n-k)));
    end
end


function [C] = control_XY(n, s)
% X and Y controls, s first qubits
    if nargin < 2
        s = n;
    end
    for k=1:s
        C{2*k-1} = mkron(speye(2^(k-1)), 0.5 * SX, speye(2^(n-k)));
        C{2*k}   = mkron(speye(2^(k-1)), 0.5 * SY, speye(2^(n-k)));
    end
end


function [C] = control_Z(n, s)
% Z controls, s first qubits
    if nargin < 2
        s = n;
    end
    for k=1:s
        C{k} = mkron(speye(2^(k-1)), 0.5 * SZ, speye(2^(n-k)));
    end
end


function H = sum_op(n, op, mult)
% Sum of operator op applied to each of n qubits.
    if nargin < 3
        mult = ones(n, 1);
    end
    
    N = 2^n;
    H = sparse(N, N);
    for k=1:n
        H = H + mkron(speye(2^(k-1)), mult(k)*op, speye(2^(n-k)));
    end
end
end

