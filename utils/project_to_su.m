function A = project_to_su(A)
% Projects A from u(n) into the traceless subalgebra, su(n).
% Essentially this eliminates the global phase.

  n = length(A);
  temp = trace(A) / n;
  A = A - temp * eye(n);
end
