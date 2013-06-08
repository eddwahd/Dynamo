function dyn = demo_resonance()
% Example: DYNAMO discovers the importance of resonant driving on its own.

% Ville Bergholm 2011-2013


%% Pauli matrices etc.

Z = [1 0; 0 -1];


%% Define the physics of the problem

% dimension vector for the quantum system
n_qubits = 1;
dim = 2 * ones(1, n_qubits);
D = prod(dim);

% Drift Hamiltonian
H_drift = -4 * Z;

% Control Hamiltonians / Liouvillians
[H_ctrl, c_labels] = control_ops(dim, 'x');

% transformed controls?
control_type = 'm';
control_par = {[-1,2]};

% pure state transfer
initial = [1 0]';
final = [0 1]';

dyn = dynamo('S ket', initial, final, H_drift, H_ctrl);
dyn.system.set_labels('Single-qubit resonant driving demo.', dim, c_labels);


%% Initial controls

% random initial controls
T = 6;
dyn.seq_init(200, T * [1, 0], control_type, control_par);
dyn.easy_control(0, 0, 1, false);


%% Now do the actual search

dyn.ui_open();
dyn.search_BFGS(dyn.full_mask(), struct('Display', 'final', 'plot_interval', 1));
dyn.analyze();
end
