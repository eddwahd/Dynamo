function dyn = fcn()
% Simulating exciton transport on a fully connected spin network.

% Ville Bergholm 2011-2012


%if nargin < 1
%    seed = floor(rand()*10000);
%end
%randseed(seed); % Optional. Allows the same pseudo-random initial controls to be generated on every run (helpful for debugging purposes, for example)


SX = [0 1; 1 0];
SY = [0 -1i; 1i 0];
SZ = [1 0; 0 -1];
I = eye(2);
SP = (SX +1i*SY)/2;


%% the physics of the problem

% Exciton transport chain with the fiducial sink as the last node,
% sucking probability from the given chain site.

n_sites = 5; % actual chain sites
dim = 2 * ones(1, n_sites+1); % dimension vector: qubits, or spin-1/2 chain
sucksite = n_sites

% 0: ground, 1: excited state, SP lowers/annihilates
a = SP;
n_op = a' * a; % == (I-SZ) / 2; % number op


%% parameters

% energy splittings
omega = [ones(1, n_sites), 0]
% site-to-site coupling
v = 1
% dephasing
gamma = zeros(1, n_sites)
% relaxation
Gamma = 1 * 1e-2 * ones(1, n_sites)
% transfer rate from sucksite to sink
suck = 1

desc = sprintf('%d-qubit FCN, suck = %g, v = %g', n_sites,suck,v);


%% cleverness

% The Hamiltonian preserves exciton number, whereas the noise
% processes we use here leave the 0-or-1 exciton manifold
% invariant. Hence we shall limit ourselves to it:

ddd = n_sites + 2; % zero- and single-exciton subspaces (sink included)
p = [1, 1+2.^(n_sites:-1:0)]; % states to keep: zero, single exciton at each site/sink
q = setdiff(1:prod(dim), p); % states to throw away

labels = {'loss', '1', '2', '3', '4', '5*', 'sink'};
% in what follows, we slice H and the Lindblad _generators_ using p
% to obtain the effective L. This works, because... TODO
% H(q,p) == 0


%% Lindblad ops

diss = cell(1, n_sites);
deph = cell(1, n_sites);
for k = 1:n_sites
    diss{k} = op_list({{sqrt(2 * Gamma(k)) * a,   k}}, dim);
    diss{k} = diss{k}(p,p);

    deph{k} = op_list({{sqrt(2 * gamma(k)) * n_op, k}}, dim);
    deph{k} = deph{k}(p,p);
end
sink = {op_list({{a, sucksite; sqrt(2 * suck) * a', n_sites+1}}, dim)};
sink{1} = sink{1}(p,p);

% Drift Liouvillian (noise / dissipation)
L_drift = superop_lindblad(sink) +superop_lindblad(diss) +superop_lindblad(deph);


%% drift Hamiltonian

J = v * 2 * [1 1 0]; % XY coupling

% fully connected graph
C = triu(ones(n_sites+1));
C(:,end) = 0; % sink not coupled to anything

H_drift = heisenberg(dim, @(s,a,b) J(s)*C(a,b))...
  +op_sum(dim, @(k) omega(k) * n_op);

H_drift = H_drift(p,p);


%% controls


% Control Hamiltonians / Liouvillians
H_ctrl = {};
for k = 1:n_sites
    temp = op_list({{sqrt(2) * n_op, k}}, dim);
    H_ctrl{end+1} = superop_lindblad({temp(p,p)}); % dephasing controls
end

%temp = op_list({{sqrt(2) * a, 1}}, dim);
%H_ctrl{3} = superop_lindblad({temp(p,p)});

% transformed controls?
control_type = char('m' + zeros(1, length(H_ctrl)));
pp = [0, 50];
control_par = {pp, pp, pp, pp, pp};



%% initial and final states

% for pure state transfer
initial = state(prod(dim(2:end)), dim); initial = initial.data; % '10..0'
final   = state(1, dim); final = final.data; % '0..01'

dyn = dynamo('SB state overlap', initial(p), final(p), H_drift, H_ctrl, L_drift);
dyn.system.set_labels(desc, labels);

% try the expensive-but-reliable gradient method
%epsilon = 1e-3;
%dyn.config.gradient_func = @(s, m) gradient_finite_diff(s, m, epsilon);


%% set up controls
T = 25;
dyn.seq_init(181, T * [0.5, 1.0], control_type, control_par);
%dyn.easy_control(1.64);
%dyn.easy_control(0);
%dyn.easy_control([4 0 0 0 2.1]);
dyn.easy_control(0 * ones(1,n_sites), 0, 0.1);
%dyn.easy_control([0.1, 0.1, 0.1], 0, 0.0, false);

%% now do the actual search

dyn.ui_open();
dyn.search_BFGS(dyn.full_mask(), struct('Display', 'final', 'plot_interval', 1));
dyn.analyze();
end
