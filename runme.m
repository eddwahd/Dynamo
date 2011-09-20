function runme(seed)
% Implementation of the different tasks from the paper
% S. Machnes et al., arXiv:1011.4874
% Ville Bergholm 2011


if nargin < 1
    seed = floor(rand()*10000);
end
randseed(seed); % Optional. Allows the same pseudo-random initial controls to be generated on every run (helpful for debugging purposes, for example)


% All definitions are in a global variable called OC
global OC; % and now we can access it too


%% Define the physics of the problem


if true
    [timeslots, T] = test_suite(20);
    squared_controls = false(1, length(OC.system.B));
else
    
n_spins = 2;
dim = 2^n_spins;


SX = [0 1; 1 0];
SY = [0 -1i; 1i 0];
SZ = [1 0; 0 -1];
I = eye(2);

SP = (SX +1i*SY)/2;

relax1 = {kron(SP, I)};
relax2 = {kron(I, SP)};
dephase1 = {kron(SZ, I)};
dephase2 = {kron(I, SZ)};
depolarize1 = {kron(SX, I), kron(SY, I), kron(SZ, I)};
depolarize2 = {kron(I, SX), kron(I, SY), kron(I, SZ)};


% Drift Hamiltonian
H_drift = (1/2) * (kron(SX,SX) + kron(SY,SY) + kron(SZ,SZ));

% Control Hamiltonians / Liouvillians
%H_ctrl = {kron(SX,I), kron(SY,I), kron(I,SX), kron(I,SY), -0.01 * superop_lindblad(dephase2)};
H_ctrl = {kron(SX,I), kron(SY,I), kron(I,SX), kron(I,SY)};

% For dissipative controls we normally want the control field to be nonnegative.
squared_controls = logical([0, 0, 0, 0]);

% noise/dissipation
% FIXME NOTE the minus sign (different convention)
L_drift = -0.002 * superop_lindblad(depolarize1);

initial = eye(dim);
final = qft(n_spins);
%final = eye(dim); 

% for task3 (pure state transfer)
%initial = [1; 0; 0; 0];
%final = [1; 0; 0; 1]/sqrt(2);


%dynamo_init('task5', initial, final, H_drift, H_ctrl, L_drift)
dynamo_init('task1', initial, final, H_drift, H_ctrl)

T = 6 * (n_spins-1) / 1; % How much time do we have to drive the system? The value specified here has empirically been shown to work well
timeslots = 20
end

%% Optimization options


if 1
    % Time-slot configuration. Can also be a vector of delta t:s.
    normal_controls = [timeslots, length(OC.system.B)];
    t_controls = [timeslots, 1];

    % Generate random initial controls
    initial_controls = 0.2 * (rand(normal_controls) - 0.5);
    if 0
        % optimize tau
        tau_par = [0.25 * T, 0.75 * T]; % min, delta

        %tau_c = 2*pi*rand(t_controls); % cosine control
        tau_c = pi/2 * ones(t_controls); % cosine control
        
        control_mask = [true(normal_controls), true(t_controls)];
    else
        % fixed tau
        tau_par = [T, 0];
        tau_c = zeros(t_controls);
        control_mask = [true(normal_controls), false(t_controls)];
    end
    initial_controls = [initial_controls, tau_c];
    
    dynamo_init_controls(initial_controls, tau_par, squared_controls);

else
    % reuse the previous, optimized controls
    %OC.seq.par = [0.25 * T, 0.75 * T];
    
    dynamo_init_controls(OC.seq.raw_controls, OC.seq.par, OC.seq.squared);

    % Which time slots do you want to modify ? All of them
    control_mask = true(size(OC.seq.raw_controls));
end


dynamo_init_opt(control_mask);


%% Now do the actual search

fprintf ('Optimizing algorithm: GRAPE (BFGS 2nd order update scheme, updating all time slices concurrently).\n\n    Please wait, this may take a while... \n\n'); drawnow;

OC.config.BFGS = struct('fminopt', struct('Display', 'off'));
termination_reason = search_BFGS();

analyze();

end

