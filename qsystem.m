classdef qsystem < matlab.mixin.Copyable
% Copyable handle class defining a quantum system.

% Ville Bergholm 2011-2012
    
  properties
    description = '' % description string
    dim              % dimension (or dim vector) of the Hilbert space of the full system S+E
    dimSE            % total dimensions of S (the part we're interested in)
                     % and the environment E, as a two-vector
    liouville        % Do the system objects X reside in Liouville or Hilbert space?
    A                % drift generator
    B                % cell vector of control generators
    B_is_Hamiltonian % true if corresponding B item corresponds to a Hamiltonian
    X_initial        % initial state
    X_final          % final state
    norm2            % squared norm of final state
    TU = []          % time unit for generators, in seconds. G = A*t/TU = A/TU * t
    state_labels   = {} % names for the Hilbert space computational basis states
    control_labels = {} % names for the controls
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
        self.A = -L_drift +1i*comm(H_drift);

        n_controls = length(H_ctrl);
        self.B = cell(1, n_controls);
        self.B_is_Hamiltonian = true(1, n_controls);

        for k=1:n_controls
            % check for Liouvillian controls
            if length(H_ctrl{k}) == self.dim
                self.B{k} = 1i*comm(H_ctrl{k}); % Hamiltonian
            else
                self.B{k} = -H_ctrl{k}; % Liouvillian
                self.B_is_Hamiltonian(k) = false;
            end
        end
    end


    function hilbert_gens(self, H_drift, H_ctrl)
    % Set up Hilbert space generators for a system.
    % (NOTE: generators are not pure Hamiltonians, there's an extra 1i!)
        self.A = 1i * H_drift;

        n_controls = length(H_ctrl);
        self.B = cell(1, n_controls);
        self.B_is_Hamiltonian = true(1, n_controls);
        for k=1:n_controls
            self.B{k} = 1i * H_ctrl{k};
        end
        %self.M = inprod_B(self.B);

        function M = inprod_B(B)
        % Computes the inner product matrix of the control operators.

            n_controls = length(B);
            M = zeros(n_controls);
            for j = 1:n_controls
                for k = 1:n_controls
                    M(j,k) = inprod(B{j}, B{k});
                end
            end
            % FIXME what about dissipative controls? superoperators?
        end
    end
  end
  
  
  methods
    function hilbert_representation(self, i, f, H_drift, H_ctrl)
    % X_* are Hilbert space objects (kets or operators)

        self.liouville = false;
        self.set_dim(i, f);
        self.X_initial = i;
        self.X_final   = f;
        self.hilbert_gens(H_drift, H_ctrl);
    end


    function vec_representation(self, i, f, H_drift, L_drift, H_ctrl)
    % X_* are Liouville space vectors corresponding to
    % vec-torized state operators.

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
        self.liouville_gens(H_drift, L_drift, H_ctrl);
    end


    function vec_gate_representation(self, i, f, H_drift, L_drift, H_ctrl)
    % X_* are Liouville space operators corresponding to
    % vec-torized unitary gates.

        self.liouville = true;
        self.set_dim(i, f);
        self.X_initial = lrmul(i, i'); % == kron(conj(i), i);
        self.X_final   = lrmul(f, f'); % == kron(conj(f), f);
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
        n_controls = length(self.B);
        
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
  end
end
