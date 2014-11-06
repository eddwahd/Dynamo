function dyn = demo_resonance_ehc(delta)
% Example: DYNAMO discovers the importance of resonant driving on its own.
%
%  dyn = demo_resonance(delta)
%
%  Optimizes a single-qubit control sequence for the state transfer |0> \to |1>.
%  delta is the detuning (difference between the carrier and
%  resonance frequencies). Uses RWA, ignores the fast-rotating term.
%
%  Resonant driving (delta = 0) yields just a constant-strength pi pulse.
%  Slightly off-resonant driving (delta = 0.1) discovers a CORPSE-type sequence.
%  Significantly off-resonant driving (delta = 5) results in a sinusoidal pulse.
    
% Ville Bergholm 2011-2014

if nargin < 1
    delta = 4;
end

%% Pauli matrices etc.

Z = [1 0; 0 -1];


%% Define the physics of the problem

% dimension vector for the quantum system
n_qubits = 1;
dim = 2 * ones(1, n_qubits);
D = prod(dim);

% Drift Hamiltonian
H_drift = -delta * Z;

% Control Hamiltonians / Liouvillians
[H_ctrl, c_labels] = control_ops(dim, 'x');

% transformed controls?
control_type = 'm';

%shi is implementing
control_par = {[-1,2]}; % u(r) = p1 + p2/2*(1-cos(r)) -- p1 = -1, p2 = 2 -- smoothness  of cosine is important
% in principle you have consecutive of these implemented
% new possibilities -- envelope functions, bandwidth limits, CRAB-type
% control

% pure state transfer
initial = [1 0]';
final = [0 1]';

dyn = dynamo('S ket', initial, final, H_drift, H_ctrl);
dyn.system.set_labels('Single-qubit resonant driving demo.', dim, c_labels);


%% Initial controls

% random initial controls
T = 6; % total time

if 1
    % start: code added by ehchen, 2014-11-05
    % prepare an ensemble of drift hamiltonians
    ensemblesize = 1;

    tempA = cell(1,ensemblesize); % a different drift generator for each ensembles? -- each element in this cell array should be a matrix... what matrix though? some hamiltonian?
    for aa = 1:length(tempA)
    %     tempA{aa} = -rand(1)*delta * Z;
        tempA{aa} = -delta * Z;
    end
    dyn.system.A = tempA;

    tempB = cell(1,ensemblesize); % control generator -- # of controls, number of ensembles -- what should be in each element of this cell array?
    for bb = 1:length(tempB);
        tempB{bb} = H_ctrl{1};
    end
    dyn.system.B = tempB;

    dyn.system.weight = ones(1,ensemblesize); % weights of each ensemble -- must be placed before seq_init (in which cache gets initialized)

    % end of code added by ehchen, 2014-11-05
end

dyn.seq_init(200, T * [1, 0], control_type, control_par);
dyn.easy_control(0, 0, 1, false);


%% Now do the actual search
dyn.ui_open();
dyn.search_BFGS(dyn.full_mask(), struct('Display', 'final', 'plot_interval', 1));
dyn.analyze();
end
