function dynamo_init(task, initial, final, H_drift, H_ctrl, L_drift)
% Initializes Dynamo for a system and an optimization task.
%
% Governing equation: \dot(X)(t) = -(A +\sum_k u_k(t) B_k) X(t) = -G(t) X(t)

% Ville Bergholm 2011

% TODO separate init of orthogonal blocks: task, system, initial/final?

% Dynamo version
version = '1.3 alpha9';

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
% Local time. TODO UTC or local time with timezone specifier would be better, but apparently MATLAB doesn't do that.
OC.config.date = datestr(now(), 31);
OC.config.task = task;

OC.config.expmFunc = @expm;

n_controls = length(H_ctrl);

% TODO temporary fix: sparse to full
H_drift = full(H_drift);
for k=1:n_controls
    H_ctrl{k} = full(H_ctrl{k});
end

input_dim = [size(initial, 2), size(final, 2)]; % check the validity of the inputs

[system, rem] = strtok(task);
[task, rem] = strtok(rem);
[phase, rem] = strtok(rem);
out = '\nOptimize a control sequence to obtain the given';
switch system
  case {'s'}
    %% Closed system S
    if nargin ~= 5
        error('Too many parameters.')
    end
    
    switch task
      case 'state'
        out = strcat(out, ' mixed state transfer');
        % TODO more efficient Hilbert space implementation?
        system_vec();
        L_drift = 0;
        system_liouville();
        OC.config.Q_func = @Q_real;
        % FIXME problems with calcPfromH_exact_gradient, eig
        % doesn't give orthonormal eigenvectors...
        
      case {'ket', 'gate'}
        if strcmp(task, 'ket')
            out = strcat(out, ' pure state transfer');
            if any(input_dim ~= 1)
                error('Initial and final states should be normalized kets.')
            end
        else
            out = strcat(out, ' unitary gate');
            if any(input_dim == 1)
                error('Initial and final states should be unitary operators.')
            end
        end

        OC.system.X_initial = initial;
        OC.system.X_final   = final;
        system_hilbert();
        
        if strcmp(phase, 'phase')
            out = strcat(out, ' (with global phase (NOTE: unphysical!))');
            OC.config.Q_func = @Q_real;
        else
            out = strcat(out, ' (ignoring global phase)');
            OC.config.Q_func = @Q_abs;
        end
        
      otherwise
        error('Unknown task.')
    end
    out = strcat(out, ' in a closed system.\n');

    % L: X_final' propagated backward 
    OC.cache.L_end = OC.system.X_final';

    % the generator is always Hermitian and thus normal => use exact gradient
    OC.config.gradientFunc = @gradient_exact;
    OC.config.calcPfromHfunc = @calcPfromH_exact_gradient; % When computing exact gradient, we get exponentiation for free due to the eigendecomposition (see paper for details)    

    
  case {'sb'}
    %% Open system S with bath B
    switch task
      case 'state'
        out = strcat(out, ' quantum state transfer');
        system_vec();

      case 'gate'
        out = strcat(out, ' quantum map');
        if any(input_dim == 1)
            error('Initial and final states should be operators.')
        end
        
        OC.system.X_initial = lrmul(initial, initial'); % == kron(conj(initial), initial);
        OC.system.X_final   = lrmul(final, final'); % == kron(conj(final), final);
        
      otherwise
        error('Unknown task.')
    end
    out = strcat(out, ' in an open system under Markovian noise.\n');
    
    system_liouville();

    % The generator isn't usually normal, so we cannot use the exact gradient method
    OC.config.Q_func = @Q_open;
    OC.opt.max_violation = 0;

    % L: reverse propagator
    OC.cache.L_end = eye(length(OC.system.X_final));
    
    %OC.config.Q_func = @Q_real; % TEST, requires also OC.cache.L{end} = X_final'
    %OC.config.gradientFunc = @gradient_first_order_aprox;
    OC.config.calcPfromHfunc = @calcPfromH_expm;

  case {'se'}
    %% Closed system S + environment E
    error('Not implemented yet.')
  
  case {'seb'}
    %% Open system S + environment E with bath B
    error('Not implemented yet.')
    
  otherwise
    error('Unknown system specification.')
end

fprintf(out);
fprintf('Optimization system dimension: %d\n', length(OC.system.X_final));

% Calculate the squared norm |X_final|^2 to scale subsequent fidelities (the real() is just taking care of rounding errors).
% We use the Hilbert-Schmidt inner product (and the induced Frobenius norm) throughout the code.
OC.system.norm2 = real(inprod(OC.system.X_final, OC.system.X_final));



function system_vec()
% Set up the vec representation for initial and final states in a Liouville space.

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


function system_hilbert()
% Set up Hilbert space generators.

  % (NOTE: generators are not pure Hamiltonians, extra 1i!)
  OC.system.A = 1i * H_drift;
  OC.system.B = cell(1, n_controls);
  OC.system.B_is_superop = false(1, n_controls);
  for k=1:n_controls
      OC.system.B{k} = 1i * H_ctrl{k};
  end
end


function system_liouville()
% Set up Liouville space generators.

  OC.system.A = L_drift +1i*comm(H_drift);
  OC.system.B = cell(1, n_controls);
  OC.system.B_is_superop = false(1, n_controls);

  % Liouville space dimension
  dim = length(OC.system.X_final);

  for k=1:n_controls
      % check for Liouvillian controls
      if length(H_ctrl{k}) ~= dim
          OC.system.B{k} = 1i*comm(H_ctrl{k}); % Hamiltonian
      else
          OC.system.B{k} = H_ctrl{k}; % Liouvillian
          OC.system.B_is_superop(k) = true;
      end
  end
end


end

