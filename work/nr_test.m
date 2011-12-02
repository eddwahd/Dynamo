function ret = nr_test(eps)
% Newton-Raphson test
% Ville Bergholm 2011
    
if nargin < 1
    eps = 1e-3;
end

V = rand_U(7); % target

U = rand_U(7);
A = logm(U);

%delta = 1i*rand_hermitian(7) * eps;
delta = (rand_hermitian(7) + 1i*rand_hermitian(7)) * eps;

[L, G] = get_LG(A);

J = U * L * (G .* (L'*delta*L)) * L';
test_J(@expm, A, delta, J);

J = L * ((L'*U'*delta*L) ./ G) * L';
test_J(@logm, U, delta, J);



% discrepancy logarithm
Q = V' * U;
F = discrepancy(U, V);

[L, G] = get_LG(F);

ret = 0;

function J = discr_J()
% Jacobian of the discrepancy
    J = L * ((L'*Q'*delta*L) ./ G) * L';
end

end


function F = discrepancy(U, V)
% fraktur L
    F = logm(V' * U);
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