function dynamo_init(task, initial, final, H_drift, H_ctrl, L_drift)
% Initializes Dynamo for a system and an optimization task.
%
% Governing equation: \dot(X)(t) = -(A +\sum_k u_k(t) B_k) X(t) = -G(t) X(t)

% Ville Bergholm 2011

% TODO separate init of orthogonal blocks: task, system, initial/final?

% Dynamo version
version = '1.3 alpha7';

% All definitions are in a global variable called OC
global OC; % and now we can access it too

if numel(OC)==0
    disp (' ');
    fprintf('DYNAMO - Quantum Dynamic Optimization Package v%s\n', version);
    disp (' ');
    disp (' ');
    disp ('(c) Shai Machnes et al. 2010-2011');
    disp ('email: shai.machnes at uni-ulm.de');
    disp (' ');
    disp ('All computer programs/code/scripts are released under the terms of the GNU Lesser General Public License 3.0 and Creative-Commons Attribution Share-Alike (see "LICENSE.txt" for details).');
    disp ('  ');
    disp ('If you use DYNAMO in your research, please add an attribution in the form of the following reference: S. Machnes et al, arXiv 1011.4874');
    disp (' ');
    disp ('For the latest version of this software, guides and information, visit http://www.qlib.info');
    disp ('  ');
    disp ('DYNAMO initialized successfully.');
    disp ('  ');    
    disp ('  ');
    drawnow;
end


task = lower(task);


%% Some basic data provenance

OC.config.version = version;
% Local time. UTC or local time with timezone specifier would be better, but apparently MATLAB doesn't do that.
OC.config.date = datestr(now(), 31);
OC.config.task = task;


n_controls = length(H_ctrl);
OC.system.B = cell(1, n_controls);
OC.system.B_is_superop = false(1, n_controls);

OC.config.expmFunc = @expm;

% We use the Hilbert-Schmidt inner product (and the induced
% Frobenius norm) throughout the code.

switch task
  case {'task1', 'task2', 'task3'}
    fprintf('Optimize a control sequence to obtain the given ')
    if strcmp(task, 'task3')
        fprintf('pure state transfer.\n')
    else
        fprintf('unitary gate.\n')
    end

    % in Hilbert space
    
    if strcmp(task, 'task2')
        fprintf('Global phase matters (unphysical, mostly a curiosity!).\n')
        OC.config.Q_func = @Q_real;
    else
        fprintf('Ignore global phase.\n')
        OC.config.Q_func = @Q_abs;
        %OC.config.Q_func = @Q_abs2;
    end
    
    % the generator is always Hermitian and thus normal => use exact gradient
    OC.config.gradientFunc = @gradient_exact;
    OC.config.calcPfromHfunc = @calcPfromH_exact_gradient; % When computing exact gradient, we get exponentiation for free due to the eigendecomposition (see paper for details)    

    % initial and target elements    
    OC.system.X_initial = initial;
    OC.system.X_final   = final;

    % (NOTE: generators are not pure Hamiltonians, extra 1i!)
    OC.system.A = 1i * H_drift;
    for k=1:n_controls
        OC.system.B{k} = 1i * H_ctrl{k};
    end

    if nargin ~= 5
        error('Too many parameters.')
    end
    
  case {'task5', 'task6'}
    fprintf ('Optimize a control sequence to obtain the given ')
    if strcmp(task, 'task5')
        fprintf('unitary quantum operation')
        OC.system.X_initial = lrmul(initial, initial'); % == kron(conj(initial), initial);
        OC.system.X_final   = lrmul(final, final'); % == kron(conj(final), final);
    else
        fprintf('quantum state transfer')
        % state vectors are converted to state operators
        if size(initial, 2) == 1
            initial = initial * initial';
        end
        if size(final, 2) == 1
            final = final * final';
        end
        
        OC.system.X_initial = vec(initial);
        OC.system.X_final   = vec(final);
    end
    fprintf(' under Markovian noise.\n');

    % Liouville space dimension
    dim = length(OC.system.X_final)
    
    % In Liouville space, vec representation (phase matters).
    % The generator isn't usually normal, so we cannot use the exact gradient method
    OC.config.Q_func = @Q_open;
    OC.opt.max_violation = 0;
    
    %OC.config.Q_func = @Q_real; % TEST, requires also OC.cache.L{end} = X_final'
    %OC.config.gradientFunc = @gradient_first_order_aprox;

    OC.config.calcPfromHfunc = @calcPfromH_expm;

    OC.system.A = L_drift +1i*comm(H_drift);
    for k=1:n_controls
        % check for Liouvillian controls
        if length(H_ctrl{k}) ~= dim
            OC.system.B{k} = 1i*comm(H_ctrl{k}); % Hamiltonian
        else
            OC.system.B{k} = H_ctrl{k}; % Liouvillian
            OC.system.B_is_superop(k) = true;
        end
    end

  otherwise
    error('Unknown task.')
    
end


% Calculate the squared norm |X_final|^2 to scale subsequent fidelities (the real() is just taking care of rounding errors)
OC.system.norm2 = real(inprod(OC.system.X_final, OC.system.X_final));

