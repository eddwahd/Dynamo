function ret = trace_matmul(A,B)
% Compute trace(A*B) efficiently.
%
% Utilizes the identity: trace(A*B) == sum(sum(transpose(A).*B)
% left side is O(n^3) to compute, right side is O(n^2)

ret = sum(sum((A.') .* B));
end
