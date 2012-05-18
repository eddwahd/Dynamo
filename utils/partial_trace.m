function res = partial_trace(x, S, q)
% Partial trace over one subsystem.
%   y = partial_trace(x, dim_S, q)
%
%   Returns y = trace_q(x), where x is an operator on H_S \otimes H_E
%   (or its vectorization). q is either 1 or 2.
%
%   y is an operator on the remaining Hilbert space.

    if S == 0
        % quick exit
        res = x;
        return
    end

    % figure out dim_E
    if size(x, 2) == 1
        % vectorized
        dim = sqrt(length(x));
    else
        % square matrix
        dim = size(x, 1);
    end
    E = dim / S;

    % build a linear index stencil for taking the partial trace
    if q == 1
        % trace over S
        R = S;
        stride = E * (dim + 1);
        stencil = zeros(E, E);
        col = 1:E;  % linear indices for the first column
        for k=0:E-1
            stencil(:, k+1) = col +(k * dim);
        end
    else
        % trace over E
        R = E;
        stride = dim + 1;
        stencil = zeros(S, S);
        col = (0:S-1) * E + 1;  % linear indices for the first column
        for k=0:S-1
            stencil(:, k+1) = col +(k * dim * E);
        end
    end

    % take the trace
    res = 0;
    for k=0:R-1
        res = res +x(stencil +k*stride);
    end

    % old code for trace_S
    %res = zeros(E, E);
    %span = 1:E;
    %for k = 1:S
    %    res = res + x(span, span);
    %    span = span + E;
    %end
end
