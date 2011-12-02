function ret = nr_test(epsilon)
% Newton-Raphson test
% Ville Bergholm 2011
    
if nargin < 1
    epsilon = 1e-3;
end

d = 4;
V = rand_U(d); % target


if false
disp('linearization tests:')
U = rand_U(d);
A = logm(U);
%delta = 1i*rand_hermitian(d) * epsilon;
delta = (rand_hermitian(d) + 1i*rand_hermitian(d)) * epsilon;

[L, G] = get_LG(A);

J = U * L * (G .* (L'*delta*L)) * L';
test_J(@expm, A, delta, J);

J = L * ((L'*U'*delta*L) ./ G) * L';
test_J(@logm, U, delta, J);
end


disp('gradient:')
S = angular_momentum(d);

f0 = randn(2,1); % parameters
delta = epsilon * randn(2,1); % parameter deltas

% discrepancy logarithm
F = discrepancy(f0);

J = build_J(f0, F);
ret = J;

cond = conditioning(vec(F), J)
Jd = J * delta;
%Jd = build_J_delta(f0, F, delta);

test_J(@discrepancy, f0, delta, inv_vec(Jd));






function F = discrepancy(f)
% fraktur L
    U = build_U(f);
    F = logm(V' * U);
end


function U = build_U(f)
% returns U corresponding to the parameters in f
  U = expm(-1i*f(2)*S{1}) * expm(-1i*f(1)*S{3});
end


function J = build_J_delta(f, F, delta)
% directional derivative (vector)

    [L, G] = get_LG(F);
    U = build_U(f);
    dU = U * (-1i*S{3}) * delta(1) +(-1i*S{1}) * U * delta(2);
    % simplification, V cancels (it still affects L and G!)
    J = vec(L * ((L'*U'* (dU) *L) ./ G) * L');
end

function J = build_J(f, F)
% full Jacobian (matrix)

    [L, G] = get_LG(F);
    U = build_U(f);
    
    for k=1:2
        if k==1
            dU = U * (-1i*S{3});
        else
            dU = (-1i*S{1}) * U;
        end
        J(:,k) = vec(L * ((L'*U'* (dU) *L) ./ G) * L');
    end
end


end % main function


function c = conditioning(L, J)
% compute the ill-conditioning
% least squares solution of L + J*p = 0
% (using the left pseudoinverse of J since J has more rows than cols)
    p = (J' * J) \ (J' * -L);
    c = norm(p);
end


function [L, G] = get_LG(A)
% Compute \Lambda and \Gamma matrices
    [L, d] = eig(A);
    d = diag(d);

    temp = kron(d, ones(size(d)).');
    temp = temp.' - temp;
    G = gggamma(temp);
end


function test_J(f, x0, delta, J)
% Try linearizing f using J, compare to exact result.
    y0 = f(x0);
    y = f(x0 + delta);
    
    norm(y-y0)
    norm(y-y0-J)
end


function ret = gggamma(z)
    ret = (exp(z) - 1) ./ z;
    ret(isnan(ret)) = 1;
end
