function Q = test(seed)
% Testing

% Ville Bergholm 2011


if nargin < 1
    seed = floor(rand()*10000);
end
randseed(seed); % Optional. Allows the same pseudo-random initial controls to be generated on every run (helpful for debugging purposes, for example)


%% Define the physics of the problem

SX = [0 1; 1 0];
SY = [0 -1i; 1i 0];
SZ = [1 0; 0 -1];
I = eye(2);
SP = (SX +1i*SY)/2;


% Exciton transport chain with the fiducial sink as the last node,
% sucking probability from the last chain site.

n_sites = 3; % actual chain sites
% dimension vector for the quantum system
dim = 2 * ones(1, n_sites+1);  % qubits, or spin-1/2 chain
D = prod(dim);

% 0: ground, 1: excited state, SP lowers/annihilates
n_op = SP' * SP; % == (I-SZ) / 2; % number op

%% parameters
omega = [1 0.5 1 0] %[rand(1, n_sites), 0];
v = 1/10;
g = [0 3 0]  %rand(1, n_sites);
G = 1e-2 * [1 1 1] %* rand(1, n_sites);
suck = 1/5 %rand();

diss = {};
deph = {};

for k = 1:n_sites
    diss = horzcat(diss, {op_list({{sqrt(G(k)) * SP,   k}}, dim)});
    deph = horzcat(deph, {op_list({{sqrt(g(k)) * n_op, k}}, dim)});
end
sink = {op_list({{SP, n_sites; sqrt(suck) * SP', n_sites+1}}, dim)};

%% FIXME HACK
N_op = op_sum(dim, n_op);
% p: permutation which sorts the excitation number in the states
[s,p] = sort(diag(N_op));


% Drift Liouvillian (noise / dissipation)
L_drift = superop_lindblad(diss) +superop_lindblad(deph) +superop_lindblad(sink);

% Drift Hamiltonian
H_drift = heisenberg(dim, @(k,s)v*(k <= n_sites-1)*(s < 3))... % chain:XY coupling
          +op_sum(dim, @(k) SZ*omega(k));

% Control Hamiltonians / Liouvillians
[H_ctrl, control_type] = control(dim, 'xy', 1);
%H_ctrl = horzcat(H_ctrl, 0.01 * superop_lindblad(dephase2));

% transformed controls?
control_par = {};
%control_par = {[], [], [0, 1]};


%initial = eye(prod(dim));
%final = qft(length(dim));
%final = eye(prod(dim));

% for pure state transfer
initial = state('1000').data;
final   = state('0001').data;


dynamo_init_control_type(control_type, control_par);
dynamo_init('SB state', initial, final, H_drift, H_ctrl, L_drift)



global OC
ddd = n_sites + 2; % zero- and single-exciton subspaces (sink included)
p = p(1:ddd); 
L = -OC.system.A;
temp = reshape(L, [D, D, D, D]);
temp = temp(p, p, p, p);
Q = reshape(temp, [ddd^2, ddd^2]);
return


%% Optimization options

control_mask = control_rand(10, 10, false);
dynamo_init_opt(control_mask);


%% Now do the actual search

fprintf('\nOptimizing algorithm: GRAPE (BFGS 2nd order update scheme, updating all time slices concurrently).\n\n    Please wait, this may take a while... \n\n'); drawnow;

% All definitions are in a global variable called OC
global OC; % and now we can access it too

termination_reason = search_BFGS(optimset('Display', 'final'));

analyze();
