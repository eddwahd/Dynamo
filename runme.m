function runme(seed)
% Simple optimization demo.

% Ville Bergholm 2011


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
dim = [2 2 2]; % two qubits
D = prod(dim);


relax1 = {kron(SP, I)};
relax2 = {kron(I, SP)};
dephase1 = {kron(SZ, I)};
dephase2 = {kron(I, SZ)};
depolarize1 = {kron(SX, I), kron(SY, I), kron(SZ, I)};
depolarize2 = {kron(I, SX), kron(I, SY), kron(I, SZ)};


% Drift Hamiltonian
H_drift = heisenberg(dim, [1 1 1]);

% Control Hamiltonians / Liouvillians
[H_ctrl, control_type] = control(dim, 'xy', 2);
%H_ctrl = horzcat(H_ctrl, 0.01 * superop_lindblad(dephase2));

% transformed controls?
control_par = {};
%control_par = {[], [], [0, 1]};


% Drift Liouvillian (noise / dissipation)
%L_drift = 0.002 * superop_lindblad(depolarize1) +0.0013 * superop_lindblad(depolarize2);

initial = eye(prod(dim));
final = qft(length(dim));
%final = eye(prod(dim));
% for pure state transfer
%initial = [1; 0; 0; 0];
%final = [1; 0; 0; 1]/sqrt(2);
% mixed states
%initial = rand_positive(D);
%final = rand_positive(D);

dynamo_init('S gate', initial, final, H_drift, H_ctrl)
%dynamo_init('SB gate', initial, final, H_drift, H_ctrl, L_drift)


dynamo_init_control_type(control_type, control_par);



%% Optimization options

control_mask = control_rand(10, 20, true, true);
dynamo_init_opt(control_mask);


%% Now do the actual search

fprintf('\nOptimizing algorithm: GRAPE (BFGS 2nd order update scheme, updating all time slices concurrently).\n\n    Please wait, this may take a while... \n\n'); drawnow;

% All definitions are in a global variable called OC
global OC; % and now we can access it too

termination_reason = search_BFGS(optimset('Display', 'final'));
%search_NR();

analyze();
