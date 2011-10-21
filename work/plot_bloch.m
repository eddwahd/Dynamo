function s = plot_bloch(linestyle, reset)
% Creates a QIT control sequence from Dynamo global data, plots the
% corresponding state evolution in a Bloch sphere.

if nargin < 2
  reset = false;

  if nargin < 1
    linestyle = 'b-';
  end
end

if reset
  figure();
  plot_bloch_sphere();
  title('State evolution under optimized sequence');
end


global OC

s.A = OC.system.A;
s.B = OC.system.B;
s.tau = OC.seq.tau;
s.control = OC.seq.control;

psi = state(inv_vec(OC.system.X_initial));

% apply sequence on state psi, plot the evolution
[out, t] = seq_propagate(psi, s, @bloch_vector);
plot_bloch_trajectory(out, linestyle);
