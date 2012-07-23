function randseed(seedvalue)
% Set the random seed.
%   randseed(seed)
%
%   A not-so-easily-repeating seed value would be for example sum(100*clock)


RandStream.setGlobalStream(RandStream('mt19937ar', 'Seed', seedvalue));
