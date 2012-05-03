function [H] = heisenberg(dim, J)
% Heisenberg spin network Hamiltonian
%
% dim is the dimension vector of the system.
% J is either a cell vector of three upper triangular coupling
% matrices, or a function handle with J(s, a, b) giving the
% coefficient of the Hamiltonian term S_sa * S_sb, where
% S_sa is the s-component of the angular momentum of spin a.
%
% Examples: J = @(s,a,b) C(s)
%
%  C = [1 1 1] gives the isotropic Heisenberg coupling.
%  C = [1 1 0] gives the XX+YY coupling.
%  C = [0 0 1] gives the Ising ZZ coupling.
    
% Ville Bergholm 2011-2012


  if isa(J, 'function_handle')
    Jfunc = J;
  else
    Jfunc = @(s, a, b) J{s}(a, b);
  end
    
  n = length(dim); % number of spins in the network
  H = sparse(0);


  for a = 1:n-1
      A = angular_momentum(dim(a)); % spin ops for first site
      for b = a+1:n
          B = angular_momentum(dim(b)); % and the second
          for s = 1:3
              H = H + op_list({{Jfunc(s, a, b) * A{s}, a; B{s}, b}}, dim);
          end
      end
  end
end
