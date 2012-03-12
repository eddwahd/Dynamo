function dyn = fmo()
% Simulating exciton transport in the FMO complex.

% Model from M.B. Plenio and S.F. Huelga, "Dephasing-assisted transport: quantum networks and biomolecules", NJP 10, 113019 (2008).
% Dipole moments etc. from F.Caruso et al, "Coherent open-loop optimal control of light-harvesting dynamics", arXiv:1103.0929v1.
    
% Ville Bergholm 2012


%randseed(825);


    
SX = [0 1; 1 0];
SY = [0 -1i; 1i 0];
SZ = [1 0; 0 -1];
I = eye(2);
SP = (SX +1i*SY)/2;


%% Define the physics of the problem

% Exciton transport network with the fiducial sink as site 8,
% sucking probability from site 3

n_sites = 7;
dim = 2 * ones(1, n_sites + 1); % dimension vector: qubits, or spin-1/2 chain

% 0: ground, 1: excited state, SP lowers/annihilates
a = SP;
n_op = a' * a; % == (I-SZ) / 2; % number op


%% parameters

% The unit for energy is 2 pi hbar c 100/m,
% and the unit for time (2 pi c 100/m)^{-1} \approx 5.31 ps.

% dephasing
gamma = 0 * ones(1, n_sites)
% optimal rates for T = 5
gamma_5 = [469.34, 5.36, 99.13, 5.55, 114.86, 1.88, 291.08];
% optimal rates for T = inf
gamma_inf = [27.4, 26.84, 1.22, 87.12, 99.59, 232.76, 88.35];
% my optimized const seq.
gamma_5_own = [22.5 158, 0, 110, 203, 162, 14];

% relaxation
Gamma = 0.5/188 * ones(1, n_sites)
% transfer rate from site 3 to sink
transfer_rate = 10/1.88

% dipole moments of the chromophores in Debye
mu = ...
  [-3.081,  2.119, -1.669;...
   -3.481, -2.083, -0.190;...
   -0.819, -3.972, -0.331;...
   -3.390,  2.111, -1.080;...
   -3.196, -2.361,  0.792;...
   -0.621,  3.636,  1.882;...
   -1.619,  2.850, -2.584];

% times maximum electric field strength in (energy unit)/D
mu = mu * 15; % 

desc = sprintf('FMO complex, transfer rate = %5.5g, Gamma = %5.5g', transfer_rate, Gamma(1));


%% cleverness

% The Hamiltonian preserves exciton number, whereas the noise
% processes we use here leave the 0-or-1 exciton manifold
% invariant. Hence we shall limit ourselves to it:
% basis: loss, 7 network sites, sink

ddd = n_sites + 2; % zero- and single-exciton subspaces (sink included)
p = [1, 1+2.^(n_sites:-1:0)]; % states to keep: zero, single exciton at each site/sink
q = setdiff(1:prod(dim), p); % states to throw away

st_labels = {'loss', '1', '2', '3*', '4', '5', '6', '7', 'sink'};


%% Lindblad ops

% NOTE the factors of sqrt(2) to match the parameters with our Lindblad eq. convention.
diss = cell(1, n_sites);
deph = cell(1, n_sites);
for k = 1:n_sites
    diss{k} = op_list({{sqrt(2 * Gamma(k)) * a,   k}}, dim);
    diss{k} = diss{k}(p,p);

    deph{k} = op_list({{sqrt(2 * gamma(k)) * n_op, k}}, dim);
    deph{k} = deph{k}(p,p);
end
sink = {op_list({{a, 3; sqrt(2 * transfer_rate) * a', n_sites+1}}, dim)};
sink{1} = sink{1}(p,p);

% Drift Liouvillian (noise / dissipation)
L_drift = superop_lindblad(sink) +superop_lindblad(diss) +superop_lindblad(deph);


%% drift Hamiltonian

H_drift = [215, -104.1,  5.1,  -4.3,   4.7, -15.1,  -7.8;
           0,      220, 32.6,   7.1,   5.4,   8.3,   0.8;
           0,        0,    0, -46.8,   1.0,  -8.1,   5.1;
           0,        0,    0,   125, -70.7, -14.7, -61.5;
           0,        0,    0,     0,   450,  89.7,  -2.5;
           0,        0,    0,     0,     0,   330,  32.7;
           0,        0,    0,     0,     0,     0,   280];
H_drift = H_drift + triu(H_drift, 1)'

% NOTE hybrid basis: rotate states 1 and 2
%U = eye(7); U(1:2,1:2) = fliplr(hadamard(2)/sqrt(2));
%H_drift = U'*H_drift*U;
H_drift = blkdiag(0, H_drift, 0);


%% controls

% Control Hamiltonians / Liouvillians
H_ctrl = {};

if 1
% dephasing controls
for k = 1:n_sites
    temp = op_list({{sqrt(2) * n_op, k}}, dim);
    H_ctrl{end+1} = superop_lindblad({temp(p,p)});
end
% transformed controls?
control_type = 'mmmmmmm';
par = [0, 500];
control_par = {par, par, par, par, par, par, par};
c_labels = {'D_1', 'D_2', 'D_3', 'D_4', 'D_5', 'D_6', 'D_7'};

else

% Laser field control ops (xr, xi, yr, yi, zr, zi)
for k = 1:3
    temp_r = 0;
    temp_i = 0;
    for s = 1:n_sites
        temp_r = temp_r -mu(s,k) * op_list({{SX, s}}, dim);
        temp_i = temp_i -mu(s,k) * op_list({{SY, s}}, dim);
    end
    H_ctrl{end+1} = temp_r(p,p); 
    H_ctrl{end+1} = temp_i(p,p); 
end
control_type = '......';
control_par = {};
c_labels = {'Xr', 'Xi', 'Yr', 'Yi', 'Zr', 'Zi'};
end


%% initial and final states

% for pure state transfer
initial = state('10000000'); initial = initial.data;
final   = state('00000001'); final = final.data;

% NOTE hybrid basis
% initial = U' * initial;
dyn = dynamo('SB state overlap', initial(p), final(p), H_drift, H_ctrl, L_drift);
dyn.system.set_labels(desc, st_labels, c_labels);


% try the expensive-but-reliable gradient method
%epsilon = 1e-3;
%dyn.config.gradient_func = @(s, m) gradient_finite_diff(s, m, epsilon);


%% set up controls

T = 5;
dyn.seq_init(109, T * [0.5, 1.0], control_type, control_par);
%dyn.easy_control(1 * gamma_5);
dyn.easy_control(0*gamma, 0, 10); % random

%return
%% now do the actual search

dyn.ui_open();
dyn.search_BFGS(dyn.full_mask(), struct('Display', 'final', 'plot_interval', 1));

dyn.analyze();
end
