function [H] = heisenberg(dim, J)
% Heisenberg spin chain Hamiltonian.
%
%  J = [1 1 1] gives the isotropic Heisenberg coupling.
%  J = [1 1 0] gives the XY coupling.
%  J = [0 0 1] gives the Ising coupling.
    
% Ville Bergholm 2011


  if isa(J, 'function_handle')
    Jfunc = J;
  else
    Jfunc = @(k,s) J(s);
  end
    
  n = length(dim); % number of spins in chain
  H = sparse(0);
  A = angular_momentum(dim(1)); % spin ops for first site

  for k=1:n-1
    % coupling between spins k and k+1: A \cdot B
    B = angular_momentum(dim(k+1)); % spin ops
    temp = Jfunc(k, 1) * kron(A{1}, B{1}) +Jfunc(k, 2) * kron(A{2}, B{2}) +Jfunc(k, 3) * kron(A{3}, B{3});
    H = H + mkron(speye(prod(dim(1:k-1))), temp, speye(prod(dim(k+2:end))));
    
    % ops for next site
    A = B;
  end
end
