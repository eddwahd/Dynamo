function dyn = demo_resonance_ehc_v4(detuning)
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

% V2 -- want to simulate ideal driving for a qutrit
% V3 -- trying to implement ensemble -- works now, but only testing with
% ensemble of size 1. it doesn't really match non-ensemble case- has to do
% with sign of weighting array?
% V4 -- trying to scale up ensemble sizes now
    
% Ville Bergholm 2011-2014
close all
if nargin < 1
    detuning = 4; % GHz
end

%% physical values
ensemblesize = 10;

we = 2*pi*0.01; % zeeman splitting, resonance of electron
we_noise = 2*pi*0.0001; % 0.1MHz amplitude noise (~0.05G B field fluctuations)
we_list = linspace(-we_noise,we_noise,ensemblesize);
rabifreq = 2*pi*0.003; % in GHz (= 10 MHz)

%% Pauli matrices etc.

Z = [1 0; 0 -1];
Z = [1 0 0; 0 0 0; 0 0 -1];
X = 1/sqrt(2)*[0 1 0; 1 0 1;0 1 0];
Ysorta = 1/sqrt(2)*[0 -1i 0; 1i 0 1i; 0 -1i 0];

weights = 

%% Define the physics of the problem

% dimension vector for the quantum system
n_qubits = 1;
dim = 3 * ones(1, n_qubits);
D = prod(dim);

% Drift Hamiltonian
H_drift = -detuning * Z;
H_drift = 2*pi*detuning * Z.^2 + we * Z;

% Control Hamiltonians / Liouvillians
% [H_ctrl, c_labels] = control_ops(dim, 'x');
H_ctrl = {X,Ysorta};
c_labels = {'x','ysorta'}; % look at derivation of RWA to see why it's not exact the Y-angular momentum matrix

% transformed controls? m means limited
control_type = 'mm';

%shi is implementing


temp = rabifreq*[-1,2]; % coefficient of control fields.. so you're allowing it to go from -rabifreq to + rabifreq
control_par = {temp, temp}; % u(r) = p1 + p2/2*(1-cos(r)) -- p1 = -1, p2 = 2 -- smoothness  of cosine is important
% in principle you have consecutive of these implemented
% new possibilities -- envelope functions, bandwidth limits, CRAB-type
% control

% pure state transfer % can 
initial = [0 1 0]';
final = [0 0 1]';

dyn = dynamo('S ket', initial, final, H_drift, H_ctrl);
% using config.error_func = @error_abs; line 168 of dynamo
% using config.gradient_func = @gradient_g_exact; line 170 of dynamo
% calculate P (propogator) from exact gradient? 'self.calcPfromHfunc = @calcP_expm_exact_gradient;'
% 'g' calculated on line 251 in cache (self.g{k} = partial_trace(temp, sys.dimSE, 1);)
dyn.system.set_labels('Single-qubit resonant driving demo.', dim, c_labels);


%% Initial controls

% random initial controls
T = 2/rabifreq; % total time

if 1
    % start: code added by ehchen, 2014-11-05
    % prepare an ensemble of drift hamiltonians

    tempA = cell(1,ensemblesize); % a different drift generator for each ensembles? -- each element in this cell array should be a matrix... what matrix though? some hamiltonian?
    for aa = 1:size(tempA,2)
    %     tempA{aa} = -rand(1)*delta * Z;
        tempA{aa} = 1i*H_drift;
%         tempA{aa} = -1i*(2*pi*detuning * Z.^2 + (we+rand(1)*we_noise) * Z);
         tempA{aa} = -1i*(2*pi*detuning * Z.^2 + (we+we_list(aa)) * Z);
    end
    dyn.system.A = tempA;

    tempB = cell(2,ensemblesize); % control generator -- # of controls, number of ensembles -- what should be in each element of this cell array?
    for bb = 1:size(tempB,2);
        tempB{1,bb} = -1i*H_ctrl{2};
        tempB{2,bb} = -1i*H_ctrl{1};
    end
    dyn.system.B = tempB;

    % weights have to be negative as opposed to not doing ensemble mode..!!!!! found the problem :) -ehchen
    dyn.system.weight = 1*ones(1,ensemblesize); % weights of each ensemble -- must be placed before seq_init (in which cache gets initialized)

    % end of code added by ehchen, 2014-11-05
end

TimeBins = 201; % should be an odd number.. usually
dyn.seq_init(TimeBins, T * [1, 0], control_type, control_par); % number of time steps, total time (end and start), 
dyn.easy_control(0, 0, 1, false); % ??


%% Now do the actual search
dyn.ui_open();
dyn.search_BFGS(dyn.full_mask(), struct('Display', 'final', 'plot_interval', 1));
dyn.analyze();
end
