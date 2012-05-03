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

dim = dyn.system.dim;

% what should we plot?
if dyn.system.liouville
    ini = state(inv_vec(dyn.system.X_initial), dim);
    fin = state(inv_vec(dyn.system.X_final), dim);
else
    ini = state(dyn.system.X_initial, dim);
    fin = state(dyn.system.X_final, dim);
end


figure()
if 1
    % apply sequence on state ini, plot the evolution
    [out, t] = seq_propagate(ini, seq, @bloch_vector, 0.01);
    plot_state_trajectory(out, linestyle, reset);

    B = {bloch_vector(fin)};
    plot_state_trajectory(B, 'm.', false);
else
    % purity plot
    [out, t] = seq_propagate(ini, seq, @purity, 0.1);
    out = real(cell2mat(out))
    dyn.plot_seq();
    hold on;
    plot(t, 2*out)
end