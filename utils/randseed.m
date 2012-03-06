function randseed(seedvalue)
% Set the random seed.

if nargin < 1
    seedvalue = sum(100*clock);
end

RandStream.setDefaultStream(RandStream('mt19937ar', 'seed', seedvalue));
end
