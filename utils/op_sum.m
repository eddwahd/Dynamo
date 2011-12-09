function H = op_sum(dim, op)
% Sum of operator op applied to each of n qubits.

if isa(op, 'function_handle')
    opfun = op;
else
    opfun = @(k) op;
end
    
n = length(dim);
G = {};
for k=1:n
    G = horzcat(G, {{opfun(k), k}});
end

H = op_list(G, dim);
