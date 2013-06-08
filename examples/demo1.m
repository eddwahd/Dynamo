function dyn = demo1(seed)
% Simple optimization demo.

% Ville Bergholm 2011-2013


if nargin < 1
    seed = floor(rand()*10000);
end
randseed(seed); % Optional. Allows the same pseudo-random initial controls to be generated on every run (helpful for debugging purposes, for example)

%% Pauli matrices

X = [0, 1; 1, 0];
Y = [0, -1i; 1i, 0];
Z = [1, 0; 0, -1];


%% Define the physics of the problem

% dimension vector for the quantum system
n_qubits = 3;
dim = 2 * ones(1, n_qubits);
D = prod(dim);

% Drift Hamiltonian
J = 2 * [0 0 1]; % Ising coupling
C = diag(ones(1, n_qubits-1), 1); % linear chain
H_drift = heisenberg(dim, @(s, a, b) J(s) * C(a,b)) -op_sum(dim, @(k) (k+2)*Z);

% Control Hamiltonians / Liouvillians
H_ctrl = {op_sum(dim, X), op_sum(dim, Y)};
c_labels = {'X', 'Y'};

% transformed controls?
control_type = '..';
control_par = {};

% Drift Liouvillian (noise / dissipation)
%L_drift = 0.002 * superop_lindblad(depolarize1) +0.0013 * superop_lindblad(depolarize2);

% gate
initial = eye(D);
final = qft(length(dim));

% pure state transfer
%initial = [1; 0; 0; 0];
%final = [1; 0; 0; 1]/sqrt(2);

% mixed state transfer
%initial = rand_positive(D);
%final = rand_positive(D);

dyn = dynamo('S gate', initial, final, H_drift, H_ctrl);
%dyn = dynamo('SB gate', initial, final, H_drift, H_ctrl, L_drift);
dyn.system.set_labels('Ising chain with non-uniform field, uniform controls.', dim, c_labels);


%% Initial controls

% random initial controls
T = 15;
dyn.seq_init(70, T * [0.5, 1.0], control_type, control_par);
dyn.easy_control([0, 0], 0, 1, true);


%% Now do the actual search

fprintf('\nUsing GRAPE (BFGS 2nd order update scheme, updating all time slices concurrently).\n\n');
dyn.ui_open();
dyn.search_BFGS(dyn.full_mask(), struct('Display', 'final', 'plot_interval', 1));
dyn.analyze();
end
