classdef qsystem < matlab.mixin.Copyable
% Copyable handle class defining a quantum system.

% Ville Bergholm 2011-2014
    
  properties
    description = ''    % description string
    dim                 % dimension (or dim vector) of the Hilbert space of the full system S+E
    dimSE               % total dimensions of S (the part we're interested in) and the environment E, as a two-vector
    liouville = false   % Do the system objects X reside in Liouville or Hilbert space?

    weight = 1          % vector, length == n_ensemble, weights for the ensemble samples
    A                   % cell vector of drift generators, size == [1, n_ensemble]
    B                   % cell array of control generators, size == [n_controls, n_ensemble]
    B_is_Hamiltonian    % vector, length == n_controls, true if corresponding B items represent Hamiltonians

    X_initial           % initial state
    X_final             % final state
    norm2               % squared norm of final state
    TU = []             % time unit for generators, in seconds. G = A*t/TU = A/TU * t
    state_labels   = {} % names for the Hilbert space computational basis states
    control_labels = {} % names for the controls
  end


  methods (Static)
      function ret = is_hermitian(H)
      % Returns true iff H is hermitian
          ret = (norm(H -H') < 1e-10);
      end
      
      function H = check_hamiltonian(H, message)
      % Raises an error if H is not a valid Hamiltonian.
          if ~qsystem.is_hermitian(H)
              error(strcat(mess, ' is not hermitian.'))
          end
          H = full(H); % eig cannot handle sparse matrices...
      end
      
      function ret = test_gens(C, n_ensemble)
      % Tests the validity of generators, converts them into
      % function handles taking the ensemble index as a parameter.

          if isa(C, 'function_handle')
              % already a function handle
              ret = C;
          else
              if iscell(C)
                  % cell array of matrices
                  s = size(C);
                  if s(1) == 1
                      ret = @(k,c) C{1,c};  % same for every ensemble member
                  elseif s(1) == n_ensemble
                      ret = @(k,c) C(k,c);
                  else
                      error('Number of rows in cell array does not match n_ensemble.')
                  end
              else
                  % single matrix (only used for A)
                  ret = @(k,c) C;
              end
          end
      end
  end


  methods (Access = private)
    function [A, B] = common_init(self, i, f, A, B)
    % Common initialization. TODO Some of this could be moved to constructor.

        % Set up the Hilbert space dimensions.
        % E may not exist, in which case it has dimension 1.
        self.dim      = size(i, 1); % S+E
        self.dimSE(1) = size(f, 1); % S
        temp = self.dim / self.dimSE(1); % E
        if floor(temp) ~= temp  % must be an integer
            error('Initial state must be an object on S+E, final state an object on S.');
        end
        self.dimSE(2) = temp;
        
        n_ensemble = self.n_ensemble();
        n_controls = size(B, 2); % FIXME assumes cell array

        % prepare the generators
        A = qsystem.test_gens(A, n_ensemble);
        B = qsystem.test_gens(B, n_ensemble);
          
        % initialize the generators
        self.A = cell(1, n_ensemble);
        self.B = cell(n_controls, n_ensemble);
        self.B_is_Hamiltonian = true(1, n_controls);
    end
  end
  
  
  methods
    function self = qsystem(weight)
    % constructor
        
        % ensemble init
        if ~isvector(weight)
            error('Ensemble weights must be given in a vector.')
        end
        self.weight = weight;
    end


    function abstract_representation(self, i, f, A, B)
    % X_ are abstract Hilbert space objects (vectors or matrices).

        [A, B] = self.common_init(i, f, A, B);
        self.X_initial = i;
        self.X_final = f;
        self.norm2 = norm2(self.X_final);

        % set the generators
        [n_controls, n_ensemble] = size(self.B);
        for k = 1:n_ensemble
            self.A{k} = A(k);
            for c = 1:n_controls
                self.B{c, k} = B(k, c);
            end
        end
    end

      
    function hilbert_representation(self, i, f, A, B, use_partial)
    % X_ are Hilbert space objects (kets or operators).
    % Used for closed system tasks: state, ket, gate, gate_partial, (TODO state_partial).
    % For _partial tasks, i \in SE, f \in S.
    % (NOTE: the generators are not pure Hamiltonians, there's an extra -1i!)
        
        [A, B] = self.common_init(i, f, A, B);
        self.X_initial = i;
        if nargin == 6 && use_partial
            % only with gate_partial
            self.X_final = kron(f, eye(self.dimSE(2)));
        else
            self.X_final = f;
        end
        
        % Calculate the squared norm |X_final|^2 to scale the fidelities with.
        % We use the Hilbert-Schmidt inner product (and the induced Frobenius norm) throughout the code.
        self.norm2 = norm2(self.X_final);
    
        % set the generators
        [n_controls, n_ensemble] = size(self.B);
        for k = 1:n_ensemble
            self.A{k} = -1i * qsystem.check_hamiltonian(A(k), sprintf('Drift Hamiltonian (ensemble %d)', k));
            for c = 1:n_controls
                self.B{c, k} = -1i * qsystem.check_hamiltonian(B(k, c), sprintf('Control Hamiltonian %d (ensemble %d)', c, k));
            end
        end
    end


    function vec_representation(self, i, f, A, B, use_states)
    % X_ are Liouville space vectors/operators corresponding to vec-torized state operators / unitary gates.
    % Used for open system tasks: state, state_partial, gate, (TODO gate_partial).

        self.liouville = true;
        [A, B] = self.common_init(i, f, A, B);
        
        if use_states
            % state vectors are converted to state operators
            if size(i, 2) == 1
                i = i * i';
            end
            if size(f, 2) == 1
                f = f * f';
            end
            self.X_initial = vec(i);
            self.X_final   = vec(f);
        else
            self.X_initial = lrmul(i, i'); % == kron(conj(i), i);
            self.X_final   = lrmul(f, f'); % == kron(conj(f), f);
        end
        self.norm2 = norm2(self.X_final);

        % Set up Liouville space generators.
        [n_controls, n_ensemble] = size(self.B)
        temp = true(n_ensemble, n_controls);
        for k = 1:n_ensemble
            self.A{k} = self.get_liouvillian(A(k), sprintf('Drift Hamiltonian (ensemble %d)', k));
            for c = 1:n_controls
                [self.B{c, k}, temp(k, c)] =...
                    self.get_liouvillian(B(k, c), sprintf('Control Hamiltonian %d (ensemble %d)', c, k));
            end
        end
        % TODO it does not make much sense to allow the type of a single
        % control to vary between ensemble members... we should maybe check it here
        self.B_is_Hamiltonian = all(temp); 
    end        


    function [C, is_H] = get_liouvillian(self, C, message)
    % Convert Hamiltonians into Liouvillians if necessary.

        if length(C) == self.dim
            % assume it's a Hamiltonian
            C = -1i * comm(qsystem.check_hamiltonian(C, message));
            is_H = true;
        else
            % assume it's a Liouvillian
            is_H = false;
        end
    end


    function set_TU(self, TU)
    % Sets the time unit for the system.
        self.TU = TU;
    end


    function set_labels(self, desc, st_labels, c_labels)
    % Describe the system, label the states and controls. The labels are cell vectors of strings.

        self.description = desc;
        D = self.dimSE(1); % label just the S states
        n_controls = size(self.B, 1);
        
        if nargin < 3 || isempty(st_labels)
            % use default state labels
            st_labels = char('0' + (1:D).');
        elseif ~iscell(st_labels)
            % it's a dim vector, use standard computational basis labeling
            dim = st_labels;
            if prod(dim) ~= prod(self.dimSE)
                error('Dimension vector given does not match the Hilbert space dimension.')
            end
            self.dim = dim; % store the dim vector (replacing the default scalar total dimension)

            % find where S ends and E starts
            n = find(cumprod(dim) == D, 1);
            ket = zeros(1,n);
            st_labels = cell(1, D);
            % build the labels
            for k=1:D
                st_labels{k} = ['|', char(ket + '0'), '\rangle'];
                for b = n:-1:1 % start from least significant digit
                    ket(b) = ket(b)+1;
                    if ket(b) < dim(b)
                        break;
                    end
                    ket(b) = 0;
                end
            end
        end

        if nargin < 4 || isempty(c_labels)
            c_labels = char('0' + (1:n_controls).');
        end
        
        if length(st_labels) ~= D
            error('Number of state labels given does not match the Hilbert space dimension of S.')
        end
        self.state_labels = st_labels;
        
        if length(c_labels) ~= n_controls
            error('Number of control labels given does not match the number of controls.')
        end
        self.control_labels = c_labels;
    end
    
    
    function ret = n_ensemble(self)
    % Returns the number of systems in the ensemble sample.
        ret = length(self.weight);
    end
  end
end
