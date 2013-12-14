function [F, p, q] = ex_to_full(H)
% Takes a Hamiltonian of an n-qubit system restricted to the
% 1-exciton manifold and converts it into the full system Hamiltonian.
% Assumes that the pairwise interactions are a combination of the forms
% (XX+YY)/2 = PM+MP  and  (XY-YX)/2 = i*(PM-MP).
% XX+YY gives a real off-diagonal element, XY-YX an imaginary one.
%
% Assumes no ZZ interactions, which could be accommodated too...

% Ville Bergholm 2012


SX = [0 1; 1 0];
SY = [0 -1i; 1i 0];
SZ = [1 0; 0 -1];
I = eye(2);
SP = (SX +1i*SY)/2;
n_op = diag([0, 1]);


n = length(H); % sites / qubits
dim = 2 * ones(1, n);


% chop up H
omega = diag(H);
R = real(triu(H, 1));
I = imag(triu(H, 1));

J = 2 * [1 1 0]; % XX+YY coupling
F = heisenberg(dim, @(s,a,b) J(s)*R(a,b))...
  +op_sum(dim, @(k) omega(k) * n_op);

% XY-YX coupling
for a = 1:n-1
    for b = a+1:n
        F = F - I(a,b)/2 * op_list({{SX, a; SY, b}, {-SY, a; SX, b}}, dim);
    end
end


p = [1+2.^(n-1:-1:0)]; % states to keep: single exciton at each site/sink
q = setdiff(1:prod(dim), p); % states to throw away
%temp = full(F(p,p))
%temp - H
end
