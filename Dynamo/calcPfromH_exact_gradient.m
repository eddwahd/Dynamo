function calcPfromH_exact_gradient(t)
% Computed once per OC.cache.H{t}, i.e. once per time slot, for a specific settings of all controls for that slot
%
%                     lambda(j) ~= lambda(k):  (exp(lambda(k)) - exp(lambda(j)))/(lambda(k)-lambda(j))
% grad_factor(j,k) =
%                     lambda(j) == lambda(k):   exp(lambda(k))
% The coefficient for <v(j) | H | v(k)> is -dT * grad_factor(j,k)
%
% References: 
%     [1] arXiv 1011.4874 (http://arxiv.org/abs/1011.4874)
%     [2] T. Levante, T. Bremi, and R. R. Ernst, J. Magn. Reson. Ser. A 121, 167 (1996)
%     [3] K. Aizu, J. Math. Phys. 4, 762 (1963)
%     [4] R. M. Wilcox, J. Math. Phys. 8, 962 (1967)

global OC;

minus_dt_H = -OC.seq.tau(t) * OC.cache.H{t};

N = length(minus_dt_H);

%% Compute the eigenvalue factors and eigenvectors of -dt*H

% TEST: trickery with 1i. Making the matrix hermitian causes eig to
% switch to an algorithm which produces orthogonal eigenvectors
% even for degenerate eigenvalues.
% TODO find an easier solution?
[v, lambda] = eig(-1i * minus_dt_H);
lambda = 1i * diag(lambda);  % note that the eigenvalues include -dt
lambdaExp = exp(lambda);

lambda_row_mat  = lambda * ones(1,N);
lambda_diff = lambda_row_mat - lambda_row_mat.'; % lambda_diff(j,k) = lambda(j) - lambda(k)

lambdaExp_row_mat = lambdaExp*ones(1,N);
lambdaExp_diff = lambdaExp_row_mat - lambdaExp_row_mat.'; % lambdaExp_diff(j,k) = exp(lambda(j)) - exp(lambda(k))

degenerate_mask = abs(lambda_diff) < 1e-10;
lambda_diff(degenerate_mask) = 1; % To prevent division by zero in next step

eig_factor = lambdaExp_diff ./ lambda_diff; % eig_factor(j,k) = (exp(lambda(j)) - exp(lambda(k)))/(lambda(j)-lambda(k))
eig_factor(degenerate_mask) = lambdaExp_row_mat(degenerate_mask); % For degenerate eigenvalues, the factor is just the exponent


OC.cache.H_v{t} = v;
OC.cache.H_eig_factor{t} = eig_factor;


%% And finally expm(-dt*H) using the eigendecomposition

OC.cache.P{t} = v * diag(lambdaExp) * v';
