function res = partial_trace(x, SE, q)
% Partial trace over one subsystem.
%   y = partial_trace(x, [dim_S dim_E], q)
%
%   x is an operator on H_S \otimes H_E (or its vectorization).
%   q is either 1 or 2, denoting which subsystem (S or E) to trace over.
%
%   y is an operator on the remaining Hilbert space
%     (or its vectorization, depending on the shape of x).

% Ville Bergholm 2012


    if SE(q) == 1
        % nothing to do, quick exit
        res = x;
        return
    end

    S = SE(1);
    E = SE(2);
    dim = prod(SE);
    if prod(size(x)) ~= dim^2
        error('The size of x does not match S*E.')
    end
    
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
            stencil(:, k+1) = col +(k * E * dim);
        end
    end
    % if x is vectorized, vectorize also the stencil.
    if size(x, 2) == 1
        stencil = stencil(:);
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
