function ret = partial_trace(x, dim1)
% Partial trace over the first subsystem.
%   y = partial_trace(x, dim_1)
%
%   Returns y = trace_1(x), where x is an operator on H_1 \otimes H_2,
%   y is an operator on H_2, and dim_1 = dim(H_1).

    dim2 = length(x) / dim1;

    span = 1:dim2;
    ret = zeros(dim2, dim2);
    for k = 1:dim1
        ret = ret + x(span, span);
        span = span + dim2;
    end
end
