function dyn = test_weights(T)
% Test ensemble optimization with weights.

% Ville Bergholm 2013



%% Pauli matrices etc.

X = [0 1; 1 0];
Y = [0 -1i; 1i 0];
Z = [1 0; 0 -1];
I = eye(2);
SP = (X +1i*Y)/2;


%% Define the physics of the problem

% dimension vector for the quantum system
n_qubits = 1;
dim = 2 * ones(1, n_qubits);
D = prod(dim);

H_drift = zeros(D);


% Control Hamiltonians / Liouvillians
[H_ctrl, c_labels] = control_ops(dim, 'xy');

% transformed controls?
control_type = 'mm';
par = [-1,2];
control_par = {par, par};

% Drift Liouvillian (noise / dissipation)
%L_drift = 0.002 * superop_lindblad(depolarize1) +0.0013 * superop_lindblad(depolarize2);

% gate
initial = eye(D);
final = R_x(pi);

% pure state transfer
%initial = [1; 0; 0; 0];
%final = [1; 0; 0; 1]/sqrt(2);

% mixed state transfer
%initial = rand_positive(D);
%final = rand_positive(D);

dyn = dynamo('S gate', initial, final, H_drift, H_ctrl);
%dyn = dynamo('SB gate', initial, final, H_drift, H_ctrl, L_drift);
dyn.system.set_labels('qubit test.', dim, c_labels);

N = 3;
tempA = dyn.system.A{1};
tempB = dyn.system.B;
n_controls = length(tempB);
dyn.system.A = cell(N,1);
dyn.system.B = cell(n_controls, N);
for k=1:N
    dyn.system.A{k} = tempA +(k-3)*0.05 * 1i*Z;
    dyn.system.B(:,k) = tempB(:);
    dyn.system.weight(k) = 1/N;
end


%% Initial controls

% random initial controls
%T = 2*pi;
dyn.seq_init(3, T * [0.1, 0.9], control_type, control_par);
dyn.easy_control(0, 0, 1, true);

tau_too = 1;

%% Now do the actual search

fprintf('\nUsing GRAPE (BFGS 2nd order update scheme, updating all time slices concurrently).\n\n    Please wait, this may take a while... \n\n'); drawnow;
dyn.ui_open();
dyn.search_BFGS(dyn.full_mask(tau_too), struct('Display', 'final', 'plot_interval', 1));

%dyn.analyze();
end
