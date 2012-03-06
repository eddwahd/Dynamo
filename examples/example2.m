% Example 2 for DPG 2012 Stuttgart

randseed(2783);

% 3-qubit Ising chain, XY controls at one end, QFT gate

q = 3;
dim = 2 * ones(1, q);

fprintf('Heisenberg chain, %d qubits, XY control at one end.\n', q)
J = 2 * [1 1 1]; % Heisenberg interaction
C = diag(ones(1, q-1), 1); % topology: linear chain
H = heisenberg(dim, @(s,a,b) J(s)*C(a,b));
C = control_ops(dim, 'xy', 1);

final = qft(q);
initial = eye(size(final));

dyn = dynamo('S gate', initial, final, H, C);
dyn.seq_init(100, 12 * [1, 0]);
dyn.easy_control(0.1 * ones(1,2));

dyn.search_BFGS(dyn.full_mask(), struct('Display', 'final', 'plot_interval', 1));
dyn.analyze();
