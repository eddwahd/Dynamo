function [C, labels] = control_ops(dim, ctrl, s)
% Returns a cell vector of control Hamiltonians (angular momentum ops).
%  [C, labels] = control_ops(dim, ctrl [, s])
%
%  dim is the system dimension vector.
%  ctrl is a string consisting of chars 'x', 'y' and 'z', denoting
%  which controls are to be generated (for all the subsystems).
%  Control operators are generated for the first s subsystems.

% Ville Bergholm 2011-2012


if nargin < 3
    s = length(dim);
end
  
ctrl = unique(lower(ctrl));
n = length(ctrl);
C = cell(s, n);
labels = {};
for k=1:s
    J = angular_momentum(dim(k));
    for j=1:n
        temp = ctrl(j) - 'x' + 1; % MATLAB indexing
        C{k,j} = mkron(speye(prod(dim(1:k-1))), J{temp}, speye(prod(dim(k+1:end))));
        labels{k, j} = sprintf('%c_%d', upper(ctrl(j)), k);
    end
end
C = C(:);
labels = labels(:);
end
