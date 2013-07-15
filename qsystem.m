classdef qsystem < matlab.mixin.Copyable
% Copyable handle class defining a quantum system.

% Ville Bergholm 2011-2013
    
  properties
    description = '' % description string
    dim              % dimension (or dim vector) of the Hilbert space of the full system S+E
    dimSE            % total dimensions of S (the part we're interested in) and the environment E, as a two-vector
    liouville        % Do the system objects X reside in Liouville or Hilbert space?

    weight = 1       % vector, length == n_ensemble, weights for the ensemble samples
    A                % cell vector of drift generators, size == [1, n_ensemble]
    B                % cell array of control generators, size == [n_controls, n_ensemble]
    B_is_Hamiltonian % vector, length == n_controls, true if corresponding B items represent Hamiltonians

    X_initial        % initial state
    X_final          % final state
    norm2            % squared norm of final state
    TU = []          % time unit for generators, in seconds. G = A*t/TU = A/TU * t
    state_labels   = {} % names for the Hilbert space computational basis states
    control_labels = {} % names for the controls
  end

  methods (Static)
      function ret = check_hamiltonian(H)
      % Returns true if H is a valid Hamiltonian.
          ret = (norm(H -H') < 1e-10);
      end
  end

  methods (Access = private)
    function set_dim(self, i, f)
    % Set up the Hilbert space dimensions.
    % E may not exist, in which case it has dimension 1.
        self.dim      = size(i, 1); % S+E
        self.dimSE(1) = size(f, 1); % S
    
        temp = self.dim / self.dimSE(1); % E
        if floor(temp) ~= temp  % must be an integer
            error('Initial state must be an object on S+E, final state an object on S.');
        end
        self.dimSE(2) = temp;
    end
      
    function liouville_gens(self, H_drift, L_drift, H_ctrl)
    % Set up Liouville space generators for a system.

        n_controls = length(H_ctrl);
        self.A = cell(1, 1);
        self.B = cell(n_controls, 1);

        self.A{1} = L_drift -1i*comm(H_drift);
        self.B_is_Hamiltonian = true(1, n_controls);

        for k=1:n_controls
            % check for Liouvillian controls
            if length(H_ctrl{k}) == self.dim
                if ~self.check_hamiltonian(H_ctrl{k})
                    error('Control Hamiltonian %d is not hermitian.', k)
                end
                self.B{k, 1} = -1i*comm(H_ctrl{k}); % Hamiltonian
            else
                self.B{k, 1} = H_ctrl{k}; % Liouvillian
                self.B_is_Hamiltonian(k) = false;
            end
        end
        self.weight = 1;
    end


    function hilbert_gens(self, H_drift, H_ctrl)
    % Set up Hilbert space generators for a system.
    % (NOTE: generators are not pure Hamiltonians, there's an extra 1i!)

        n_controls = length(H_ctrl);
        self.A = cell(1, 1);
        self.B = cell(n_controls, 1);
        
        if ~self.check_hamiltonian(H_drift)
            error('The drift Hamiltonian is not hermitian.')
        end
        self.A{1} = -1i * H_drift;
        self.B_is_Hamiltonian = true(1, n_controls);
        for k=1:n_controls
            if ~self.check_hamiltonian(H_ctrl{k})
                error('Control Hamiltonian %d is not hermitian.', k)
            end
            self.B{k, 1} = -1i * H_ctrl{k};
        end
        self.weight = 1;
        %self.M = inprod_B(self.B);

        function M = inprod_B(B)
        % Computes the inner product matrix of the control operators.

            M = zeros(n_controls);
            for j = 1:n_controls
                for k = 1:n_controls
                    M(j,k) = inprod(B{j}, B{k});
                end
            end
            % FIXME what about dissipative controls? superoperators?
        end
    end
  
    function abstract_gens(self, A, B)
    % Set up abstract vector space generators for a system.

        n_controls = length(B);
        self.A = cell(1, 1);
        self.B = cell(n_controls, 1);
        
        self.A{1} = A;
        self.B_is_Hamiltonian = true(1, n_controls);
        for k=1:n_controls
            self.B{k, 1} = B{k};
        end
        self.weight = 1;
    end
  end
  
  
  methods
    function abstract_representation(self, i, f, A, B)
    % X_ are Hilbert space objects (vectors or matrices).
        
        self.liouville = false;
        self.set_dim(i, f);
        self.X_initial = i;
        self.X_final = f;
        
        % Calculate the squared norm |X_final|^2 to scale the fidelities with.
        % We use the Hilbert-Schmidt inner product (and the induced Frobenius norm) throughout the code.
        self.norm2 = norm2(self.X_final);
        self.abstract_gens(A, B);
    end

      
    function hilbert_representation(self, i, f, H_drift, H_ctrl, enlarge_f)
    % X_ are Hilbert space objects (kets or operators).
    % closed system: state, ket, gate, gate_partial, (TODO state_partial)
    % For _partial tasks, i \in SE, f \in S.
        
        self.liouville = false;
        self.set_dim(i, f);
        self.X_initial = i;
        if nargin == 6 && enlarge_f
            % only with gate_partial
            self.X_final = kron(f, eye(self.dimSE(2)));
        else
            self.X_final = f;
        end
        
        % Calculate the squared norm |X_final|^2 to scale the fidelities with.
        % We use the Hilbert-Schmidt inner product (and the induced Frobenius norm) throughout the code.
        self.norm2 = norm2(self.X_final);
        self.hilbert_gens(H_drift, H_ctrl);
    end


    function vec_representation(self, i, f, H_drift, L_drift, H_ctrl)
    % X_* are Liouville space vectors corresponding to vec-torized state operators.
    % open system: state, state_partial
        
        self.liouville = true;
        self.set_dim(i, f);
        % state vectors are converted to state operators
        if size(i, 2) == 1
            i = i * i';
        end
        if size(f, 2) == 1
            f = f * f';
        end
        self.X_initial = vec(i);
        self.X_final   = vec(f);
        self.norm2 = norm2(self.X_final);
        self.liouville_gens(H_drift, L_drift, H_ctrl);
    end


    function vec_gate_representation(self, i, f, H_drift, L_drift, H_ctrl)
    % X_* are Liouville space operators corresponding to vec-torized unitary gates.
    % open system: gate
        
        self.liouville = true;
        self.set_dim(i, f);
        self.X_initial = lrmul(i, i'); % == kron(conj(i), i);
        self.X_final   = lrmul(f, f'); % == kron(conj(f), f);
        self.norm2 = norm2(self.X_final);
        self.liouville_gens(H_drift, L_drift, H_ctrl);
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
