function dyn = runme(seed)
% Simple optimization demo.

% Ville Bergholm 2011-2012


if nargin < 1
    seed = floor(rand()*10000);
end
randseed(seed); % Optional. Allows the same pseudo-random initial controls to be generated on every run (helpful for debugging purposes, for example)


%% Pauli matrices etc.

SX = [0 1; 1 0];
SY = [0 -1i; 1i 0];
SZ = [1 0; 0 -1];
I = eye(2);
SP = (SX +1i*SY)/2;


%% Define the physics of the problem

% dimension vector for the quantum system
dim = [2 2 2]; % n qubits
D = prod(dim);


relax1 = {kron(SP, I)};
relax2 = {kron(I, SP)};
dephase1 = {kron(SZ, I)};
dephase2 = {kron(I, SZ)};
depolarize1 = {kron(SX, I), kron(SY, I), kron(SZ, I)};
depolarize2 = {kron(I, SX), kron(I, SY), kron(I, SZ)};


% Drift Hamiltonian
H_drift = heisenberg(dim, [0 0 4]) -op_sum(dim, @(k) (k+2)*SZ);

% Control Hamiltonians / Liouvillians
%[H_ctrl, control_type] = control(dim, 'xy', 2);
H_ctrl = {op_sum(dim, @(k) SX), op_sum(dim, @(k) SY)};
control_type = '..';

% transformed controls?
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


%% Initial controls

% random initial controls
control_mask = dyn.init_control(125, 100, false, [], control_type, control_par);


%% Now do the actual search

fprintf('\nUsing GRAPE (BFGS 2nd order update scheme, updating all time slices concurrently).\n\n    Please wait, this may take a while... \n\n'); drawnow;
termination_reason = dyn.search_BFGS(control_mask,...
    struct('Display', 'final', 'plot_interval', 1));
dyn.analyze();
end
