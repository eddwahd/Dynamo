function seq = plot_bloch(dyn, linestyle, reset)
% Creates a QIT control sequence from Dynamo data, plots the
% corresponding state evolution in a Bloch sphere.

if nargin < 3
  reset = true;

  if nargin < 2
    linestyle = 'b-';
  end
end


% extract the sequence
seq.A = dyn.system.A;
seq.B = dyn.system.B;
seq.tau = dyn.seq.tau;
seq.control = dyn.seq.fields;

dim = [2 2];

psi = state(inv_vec(dyn.system.X_initial), dim);

figure()
if 1
    % apply sequence on state psi, plot the evolution
    [out, t] = seq_propagate(psi, seq, @bloch_vector, 0.01);
    plot_state_trajectory(out, linestyle, reset);

    s = state(inv_vec(dyn.system.X_final), dim);
    B = {bloch_vector(s)};
    plot_state_trajectory(B, 'm.', false);
else
    % purity plot
    [out, t] = seq_propagate(psi, seq, @purity, 0.1);
    out = real(cell2mat(out))
    dyn.plot_seq();
    hold on;
    plot(t, 2*out)
end