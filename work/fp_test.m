function [L, A, H_break] = fp_test(kappa, phi, tol, H)
% Construct a Lindbladian, find its fixed point states.

% Ville Bergholm 2011-2012

    
if nargin < 3
    tol = 1e-2;
end
if nargin < 4
    H = 0;
end
do_plot = true;


[L, scale] = lindblad(kappa, phi, H);
return;
[A, spectrum] = lind_fp(L, tol);

svd_tol = 5e-2;

if do_plot
    % plot a random pure state trajectory

    %s = u_propagate(state('00', [2 2]), kron(rand_U(2), rand_U(2)));
    s = u_propagate(state('00', [2 2]), rand_U(4));
    %s = state('00');
    %s = normalize(state('00')+state('10'));
    
    N = 300;
    t = linspace(0, 30/scale, N);
    out = propagate(s, L, t, @(s,H) bloch_vector(s));
    figure()
    plot_state_trajectory2(out, 'b-');
end

%F = fidelity(sFP, bloch_state(out{end}, [2 2]))

% plot the fixed point set
switch length(A)
  case 1
    disp('unique FP')
    s = state(A{1}, [2 2])
    C = concurrence(s)
    B = {bloch_vector(s)};
    
    % find a Hamiltonian which does not commute with s
    S = nullspace_hermitian(comm(s.data), svd_tol);
    [th, u, v] = principal_angles(eye(size(S, 1)), S);
    th
    H = inv_vec(u(:,end));
    
  case 2
    disp('1D FP set')
    q = sqrt(1 - trace(A{1}^2));

    s1 = state(A{1} -q*A{2}, [2 2])
    eig(s1.data)
    C1 = concurrence(s1)

    s2 = state(A{1} +q/3*A{2}, [2 2])
    eig(s2.data)
    C2 = concurrence(s2)

    B = {bloch_vector(s1), bloch_vector(s2)};

    % find the null space of [rho, *] within the real vector space of hermitian matrices
    S1 = nullspace_hermitian(comm(s1.data), svd_tol);
    S2 = nullspace_hermitian(comm(s2.data), svd_tol);
    %I = subspace_intersection(S1, S2, svd_tol);
    fprintf('Dims of commuting Hamiltonian subspaces: s1: %d, s2: %d\n', size(S1, 2), size(S2, 2));
    [th, u, v] = principal_angles(S1, S2);
    th
    % Hamiltonian which commutes with s2 but not with s1...
    H = inv_vec(u(:,end));

    
  otherwise
    error('too many FPs')
end

H_break = H + H'; % fix num. errors
if do_plot
    plot_state_trajectory2(B, 'r:o', false);
end
end



function [L, scale] = lindblad(kappa, phi, H_add)
% Construct a Lindbladian.

X = [0 1; 1 0];
Y = [0 -1i; 1i 0];
Z = [1 0; 0 -1];
P1 = [0 0; 0 1];
I = eye(2);
SP = (X +1i*Y)/2;

a = SP;     % annihilation op: |0><1|
n = a' * a; % number op


if nargin < 1
  kappa = 5;
end

H = 0;
scale = 1;

switch 'vmc_distillation'
  case 'muschik'
    % Christine Muschik's system
    % FP at |kappa| \to \infty: |00> -sgn(kappa)*|11>
    % if |kappa| >> 0, another quasi-FP appears in the middle of the
    % opposing tetrahedron face.

    % annihilation op for each qubit
    p = kron(a, I);
    q = kron(I, a);

    mu = cosh(kappa);
    %p = exp(2i * pi * rand()); % relative phase
    nu = exp(1i * phi) * sinh(kappa);
    A = {mu*p + nu*q', mu*q + nu*p'};

    scale = mu^2;
    H = kron(Z, Z);
    %H = zeros(4);
    %H = kron(X, X) +kron(Y, Y) +kron(Z, Z);
    %H = kron(X, X) +kron(Y, Y);
    %H = kron(Z, Z);
    %H = rand_hermitian(4);

  case 'muschik_variant'
    % normal instead of hyperbolic functions
    
    % annihilation op for each qubit
    p = kron(a, I);
    q = kron(I, a);
    
    mu = cos(kappa);
    nu = exp(1i * phi) * sin(kappa);

    F = mu * p +nu * q';
    G = mu * q +nu * p';
    A = {F, G};
    %H = kron(Z, Z);
    
  case 'vmc_distillation'
    % Vollbrecht-Muschik-Cirac continuous entanglement distillation procedure
    % Six qubits divided in three pairs, each pair shared between Alice and Bob.
    dim = 2 * ones(1, 6); % order: A:T,s1,s2; B:T,s1,s2 
    D = prod(dim);
    
    kappa = 0.77;
    mu = cosh(kappa);
    nu = exp(1i * phi) * sinh(kappa);
    gamma = 1;
    e_c = 0.0001;
    e_h = 0.0001;
    e_d = 0.0001;
    delta = 0.1;

    L = sparse(D^2, D^2);
    
    % each source qubit pair:
    for k = 2:3
        A_ent = {op_list({{mu * a, k}, {nu * a', 3+k}}, dim),...
                 op_list({{mu * a, 3+k}, {nu * a', k}}, dim)}; % entangling

        A_cool = {op_list({{a,  k}}, dim),   op_list({{a,  3+k}}, dim)}; % cooling
        A_heat = {op_list({{a', k}}, dim),   op_list({{a', 3+k}}, dim)}; % heating
        A_deph = {op_list({{n,  k}}, dim),   op_list({{n,  3+k}}, dim)}; % dephasing

        L = L +gamma * superop_lindblad(A_ent) +e_c * superop_lindblad(A_cool)...
            +e_h * superop_lindblad(A_heat) +e_d * superop_lindblad(A_deph);
    end 

    % annihilation op for each qubit in a pair:
    %p = kron(a, I);
    %q = kron(I, a);
    %    A = {mu * p + nu * q', mu * q + nu * p',... % entangling
    %     e_c * p,          e_c * q, ...         % cooling
    %     e_h * p',         e_h * q',...         % heating
    %     e_d * kron(n, I), e_d * kron(I, n)};   % dephasing
    %L_s = lmap(superop_lindblad(A), {[4 4], [4 4]});
    %L = gate.two(L_s, [2, 5], dim) +gate.two(L_s, [3, 6], dim);

    % distilling Hamiltonian
    s00 = sparse(state('001', [2 2 2]).data);
    s01 = sparse(state('010', [2 2 2]).data);
    s10 = sparse(state('101', [2 2 2]).data);
    s11 = sparse(state('110', [2 2 2]).data);
    F = s00*s00' +s10*s01' +s01*s10' +s11*s11';
    FA = kron(F, speye(8));
    FB = kron(speye(8), F);
    
    L = L +superop_lindblad({}, delta * (-FA +FB));
    return
    
  case 'fu'
    % Fu-Li-Li NV/resonator system

    g = kappa;
    scale = g;
    
    x = sqrt(g) / 5;
    A = {x * (p' + q')};
    Delta = g / 100;
    Theta = -g * sqrt(2)/10;
    %H = Delta*(kron(P1, I)-kron(I, P1)) -Theta*(kron(X, I) +kron(I, X)); % Li eq (11)
    H = 2*Delta*(kron(P1, I)-kron(I, P1)) +Theta*(kron(X, I) +kron(I, X)); % Li eq (15)

    sFP = normalize(Delta * state('11') +Theta * sqrt(2) * state('bell1'));
    
  case 'singlet'
    % pure singlet as FP
    x = 1/sqrt(2);
    U = [0, 1, 0, 0; x, 0, x, 0; -x, 0, x, 0; 0 0 0 1];
    %U = rand_U(4);
    A = {U*p*U', U*q*U'};

  case '1q'
    % single qubit
    A = {mu*SP + nu*SP'}
    %A = {SP};
    %H = 2*X;

  otherwise
    error('Unknown system.')
end

L = superop_lindblad(A, scale * (H + H_add));
%figure(); spy(L)
end



function rho = find_fp(L)
% hermitian fixed point of L

    [U,S,V] = svds(L, 8, 0.00001);
    diag(S)
    r = inv_vec(V(:, end));
    rho = (r + r')/2;
    rho = rho/trace(rho);
    norm(L'*L*vec(rho))
    s = state(rho, 2*ones(1,6));
    t = ptrace(s, [2 3 5 6]);
    concurrence(t)
    s1 = ptrace(s, [1 3 4 6]);
    concurrence(s1)
end