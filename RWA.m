function [Q, Q_slow, a_min, a_max, omega_rot] = RWA(H0, A, carrier, verbosity, allow_crosstalk, safety_factor)
% RWA  Rotating wave approximation.
%  [Q_static, Q_slow, a_min, a_max, omega_rot] = RWA(H0, A [, omega_carrier])
%
%  Performs the rotating frame transformation U = \exp(i H0 t), where
%  H0 is the hermitian generator of the transformation, on the
%  hermitian operator A, giving A' = U A U^\dagger.
%  The fast-rotating terms are then thrown away.
%
%  The matrix Q_static is a sum of the (nearly) static terms of A'.
%
%  Q_slow is a cell vector of function handles @(t) representing the
%  terms rotating with a slow but nonvanishing rate.
%
%  a_min is a vector, each element denoting the
%  control amplitude above which the corresponding term in Q_slow
%  first passes the safety_factor limit and thus becomes significant.
%  The higher the a_min, the less significant the term.
%
%  If omega_carrier is given, A is assumed to oscillate in time,
%
%    A(t) = A a(t) \cos(\omega_carrier t +\phi(t)). 
%
%  The maximum amplitude a_max >= |a(t)| is chosen low enough so that only the
%  transition(s) closest to omega_carrier are excited (no crosstalk), 
%  and the fast-rotating terms can be neglected.
%  In this case A' = a(t) Q(t, phi).
%  Q_static is a function handle @(phi) and Q_slow a cell
%  vector of handles @(t, phi).
%
%  The inputs are dimensionless (Hamiltonians divided by \hbar and multiplied by TU,
%  angular frequencies multiplied by TU).

% Ville Bergholm 2013-2014


% The logic here is the following:
%
% For each H0 eigenspace there is a corresponding eigenfrequency omega_k.
% Each (a,b) "block" of A rotates with the frequency delta_ab = omega_a -omega_b.
% cos(omega_carrier t) consists of two counter-rotating terms.
% Combined with delta_ab this yields a faster and a slower
% component, |omega_fast| >= |omega_slow| (equal if omega_carrier = 0 or a = b).
% The term exp(i omega_xxxx t) a_max A_ab can be neglected if the total rotation speed
% |omega_xxxx| > safety_factor * |A_ab| * a_max.
% To obtain a fixed A', we must be able to ignore all the faster components,
% and the slower components of unwanted transitions (crosstalk), which is
% done by setting a_max low enough.


% omega-tolerance for combining rotating terms/blocks/modes
tol_combine = 1e-3;
% |A|-tolerance for keeping block groups
tol_zeronorm = 0;
% relative omega-tolerance for deciding which modes are unwanted
% and, if crosstalk is to be avoided, included in a_max determination
tol_target = 1e-7;
% ignored terms with |omega|/|A| ratios above this will not be even shown
tol_show = 1e5;
% |omega|-tolerance for considering terms static
tol_veryslow = 1e-3;


format = ' %s:  2pi*|%6.4g| / %6.4g = %6.4g\n';


if nargin < 6
    safety_factor = 500; % TODO justify
    
if nargin < 5
    allow_crosstalk = false;

if nargin < 4
    verbosity = 1;

% If carrier > 0 is given, H is multiplied by cos(carrier*t).
if nargin < 3
    carrier = 0;
elseif carrier <= 0
    error('Carrier frequency must be positive.')
end
end
end
end

% The identity component of a Hamiltonian A plays no physical role,
% but affects the op_norm calculation later on when a == b.
% Removing the id. comp. always lowers the Frobenius norm.
dim = length(A);
trace_A = trace(A);
A = A -trace_A * eye(dim) / dim;

% ascending eigenenergies
[omega, P] = spectral_decomposition(H0, true);
n = length(omega);


%% first pass: compute the slow and fast rotation freqs for each block

M = n*(n+1)/2;
freqs = NaN(M, 4);

% loop over the diagonal and upper triangle of eigenspaces
k = 1;
for a=1:n
    for b=a:n
        % frame rotation angular frequency
        delta = omega(a) -omega(b); % always < 0 since omegas are sorted
        
        % slower rotating mode
        omega_slow = carrier +delta;
        % faster rotating mode
        omega_fast = carrier -delta;
        
        freqs(k, :) = [a, b, omega_slow, omega_fast];
        k = k+1;
    end
end

% sort the blocks according to omega_slow
% (if two omega_slow:s are equal, so are the corresponding omega_fast:s)
[~, ind] = sort(freqs(:, 3));
freqs = freqs(ind, :);


%% second pass: combine together all blocks with the same omega_slow/omega_fast
% All omega_slow <= carrier <= omega_fast, with equality only for diagonal blocks.
% Hence by grouping by omega_slow, we will never mix diagonal and nondiagonal blocks.

groups = zeros(0, 6);
group_label = {};
group_proj = {};
j_slow = 1;
j_fast = 2;
j_A_norm = 3;
j_vol_slow = 4;
j_vol_fast = 5;
j_diagonal = 6;
skipped = '';  % zero-norm groups skipped

if verbosity >= 2
    figure();
end
    
temp = 0;
label = '';
A_proj = 0;  % sum of block projectors
omega_slow = freqs(1, 3);
for k=1:M
    a = freqs(k, 1);
    b = freqs(k, 2);
    A_proj = A_proj +P{a} * A * P{b};
    label = sprintf('%s (%d,%d)', label, a, b);
    
    if k+1 <= M
        % still more rows to process
        temp = freqs(k+1, 3);
        if abs(temp -omega_slow) < tol_combine
            % same omega_slow, keep summing
            continue
        end
    end
    % omega_slow changes value, or we ran out of rows
    
    % norm of the summed A projection
    A_norm = op_norm(A_proj);

    omega_fast = freqs(k, 4);
    
    % TODO with carrier > 0, drop diagonal blocks since they will
    % be discarded later anyway?
    
    if A_norm > tol_zeronorm
        % ignore terms with a tiny A_norm, they have no effect on the result    

        if a ~= b
            % the "lower triangle" (+h.c.) terms counts towards the norm, too
            A_norm = A_norm * sqrt(2);
            if carrier ~= 0
                % slow or fast mode alone (for carrier == 0 they're the same)
                A_norm = A_norm * 0.5;
            end
        end
            
        % record the block group
        groups(end+1, :) = [omega_slow, omega_fast, A_norm, abs(omega_slow)/A_norm, abs(omega_fast)/A_norm, a == b];
        group_label{end+1} = label;
        group_proj{end+1} = A_proj;
        
        if verbosity >= 2
            semilogy(omega_fast/2/pi, A_norm, 'r+', abs(omega_slow)/2/pi, A_norm, 'bo')
            hold on
            text(abs(omega_slow)/2/pi, A_norm, label);
        end
    else
        % this group does nothing, just make a note of it
        skipped = [skipped, label];
    end
    
    % next group
    label = '';
    A_proj = 0;
    omega_slow = temp;
end


%% in the carrier mode, determine a_max such that the fast mode(s) (and maybe crosstalk) do not interfere

if carrier ~= 0
    % eliminate fast modes
    % NOTE On the diagonal |omega_slow| == |omega_fast|, so this
    % also eliminates diagonal blocks entirely when there is a carrier.
    % We thus cannot target the diagonal blocks... which would be crazy anyway?
    [vol_fast, j] = min(groups(:, j_vol_fast));
    
    % find lowest |omega_slow|
    abs_omega_slow = abs(groups(:, j_slow));
    temp = min(abs_omega_slow);

    % all the slower modes within tol_target*carrier of the lowest one are the 'targeted' ones
    targets = find(abs_omega_slow <= temp +tol_target * carrier);
    unwanted = setdiff(1:length(abs_omega_slow), targets);
    
    if allow_crosstalk
        vol_slow = inf; k = 1; % do not care about crosstalk
    else
        % eliminate crosstalk: all slow modes except the targeted
        % ones must have negligible effect.
        [vol_slow, k] = min(groups(unwanted, j_vol_slow));
    end

    if verbosity >= 1
        fprintf('\nf_carrier: %g, targets: ', carrier/2/pi);
        fprintf('%s, ', group_label{targets});
        fprintf('\n');
        if ~allow_crosstalk
            fprintf('crosstalk limit %s: %g\n', group_label{unwanted(k)}, vol_slow/safety_factor);
        end
	% safety factor FIXME not printed!!
        fprintf(['fast mode limit', format], group_label{j}, groups(j, j_slow)/2/pi, groups(j, j_A_norm), vol_fast/safety_factor);
    end
    a_max = min(vol_fast, vol_slow) / safety_factor;
    Q = @(phi) 0;
else
    a_max = 1; % there is no amplitude
    if verbosity >= 1
        fprintf('\nf_carrier: %g\n', carrier/2/pi);
    end
    Q = 0;
end

if verbosity >= 2
    grid on
    a = axis();
    line([1, 1] * carrier/2/pi, [a(3), a(4)]);
    temp = linspace(1, a(2), 500);
    line(temp, temp * 2*pi/(safety_factor*a_max));
    xlabel('f (1/TU)')
    ylabel('|P_a A P_b|')
    title(sprintf('f_{carrier}: %g, a_{max}: %g', carrier/2/pi, a_max))
end


%% third pass: see which slowly rotating terms to keep

Q_slow = {};
a_min = [];
omega_rot = [];

noshow = '';

% loop over the groups, in order of ascending vol_slow
[~, ind] = sort(groups(:, j_vol_slow));
for k=ind.'
    label = group_label{k};
    omega_slow = groups(k, j_slow);
    % NOTE the a_max: assume we are driving at max amplitude
    A_norm = groups(k, j_A_norm) * a_max;
    % "volatility" of the term. The higher, the less important.
    vol = abs(omega_slow) / A_norm;

    f_slow = omega_slow / (2*pi);
    
    if vol >= 0.9999 * safety_factor; % numerical safety
        % Even the slower mode is too fast, ignore.
        if verbosity >= 1
            if vol <= tol_show
                fprintf(['ignore  ', format], label, f_slow, A_norm, vol);
            else
                noshow = [noshow, label];
            end
        end
        continue
    end
    % omega_slow term is slow enough to keep.

    % build the rotating term
    A_proj = group_proj{k};
    if carrier == 0
        % slow and fast modes are the same, no phi
        R = @(t) exp(1i*omega_slow*t) * A_proj;
    else
        % ignore the fast mode (a_max takes care of it), hence the 0.5
        % diagonal modes have also been discarded
        R = @(t, phi) 0.5 * exp(1i*(omega_slow*t +phi)) * A_proj;
    end

    % In pass one we only looped over the "upper triangle" of projectors, so
    if ~groups(k, j_diagonal)
        % off-diagonal group, add the hermitian conjugate
        S = @(x) x+x';
    else
        % do nothing
        S = @(x) x;
    end

    if abs(omega_slow) <= tol_veryslow
        % Very slow component wrt. TU, add directly to Q.
        if verbosity >= 1
            fprintf(['fixed   ', format], label, f_slow, A_norm, vol);
        end
        if carrier == 0
            Q = Q +S(R(0));
        else
            % HACK recursive function definition... evaluating
            % it and saving the results before use improves performance.
            Q = @(phi) Q(phi) +S(R(0, phi));
        end
    else
        % Slowly rotating component.
        % Can happen anywhere with carrier > 0, or with carrier = 0 in off-diagonal groups.
            
        % amplitude at which the volatility of the term exceeds the safety_factor
        a_min(end+1) = groups(k, j_vol_slow) / safety_factor;
        % rotation speed
        omega_rot(end+1) = omega_slow;
        if verbosity >= 1
            fprintf(['rotating', format], label, f_slow, A_norm, vol);
        end
        if carrier == 0
            Q_slow{end+1} = @(t) S(R(t));
        else
            Q_slow{end+1} = @(t, phi) S(R(t, phi));
        end
    end
end

if ~isempty(noshow)
    fprintf('ignore+  %s\n', noshow);
end
if ~isempty(skipped)
    fprintf('zero norm%s\n', skipped);
end
end


function ret = op_norm(A)
% Given a traceless Hermitian operator A, returns a positive scalar describing the
% magnitude of the unitary flow generated by it.
  ret = norm(A, 'fro');
end
