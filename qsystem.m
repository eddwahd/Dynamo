classdef qsystem < matlab.mixin.Copyable
% Copyable handle class defining a quantum system.

% Ville Bergholm 2011-2012
    
  properties
    description = '' % System description string
    liouvillian      % Are the objects below in Liouville space or Hilbert space?
    A                % drift generator
    B                % cell vector of control generators
    B_is_Hamiltonian % true if corresponding B item is a Hamiltonian
    M
    X_initial        % initial state
    X_final          % final state
    norm2            % squared norm of final state
    state_labels   = {} % names for the Hilbert space computational basis states
    control_labels = {} % names for the controls
  end
    
  methods
    function std_representation(self, i, f)
    % X_* are Hilbert space objects (kets or operators)
        self.X_initial = i;
        self.X_final   = f;
    end


    function vec_representation(self, i, f)
    % X_* are Liouville space vectors corresponding to
    % vec-torized state operators.

    % state vectors are converted to state operators
        if size(i, 2) == 1
            i = i * i';
        end
        if size(f, 2) == 1
            f = f * f';
        end

        self.X_initial = vec(i);
        self.X_final   = vec(f);
    end


    function vec_gate_representation(self, i, f)
    % X_* are Liouville space operators corresponding to
    % vec-torized unitary gates.
        self.X_initial = lrmul(i, i'); % == kron(conj(i), i);
        self.X_final   = lrmul(f, f'); % == kron(conj(f), f);
    end


    function hilbert(self, H_drift, H_ctrl)
    % Set up Hilbert space generators for a system.

        self.liouvillian = false;
        % (NOTE: generators are not pure Hamiltonians, there's an extra 1i!)
        self.A = 1i * H_drift;

        n_controls = length(H_ctrl);
        self.B = cell(1, n_controls);
        self.B_is_Hamiltonian = true(1, n_controls);
        for k=1:n_controls
            self.B{k} = 1i * H_ctrl{k};
        end
        self.M = inprod_B(self.B);

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


    function liouville(self, H_drift, L_drift, H_ctrl)
    % Set up Liouville space generators for a system.

        self.liouvillian = true;
        self.A = -L_drift +1i*comm(H_drift);

        n_controls = length(H_ctrl);
        self.B = cell(1, n_controls);
        self.B_is_Hamiltonian = true(1, n_controls);

        % Liouville space dimension
        dim = length(self.X_final);

        for k=1:n_controls
            % check for Liouvillian controls
            if length(H_ctrl{k}) ~= dim
                self.B{k} = 1i*comm(H_ctrl{k}); % Hamiltonian
            else
                self.B{k} = -H_ctrl{k}; % Liouvillian
                self.B_is_Hamiltonian(k) = false;
            end
        end
        self.M = eye(n_controls); % FIXME temporary fix, meaningless
    end


    function ret = dimension(self)
    % Returns the dimension of the Hilbert space of the system.

        temp = length(self.A);
        if self.liouvillian
            ret = sqrt(temp);
        else
            ret = temp;
        end
    end


    function set_labels(self, desc, st_labels, c_labels)
    % Describe the system, label the states and controls. The labels are cell vectors of strings.

        self.description = desc;
        
        dim = self.dimension();
        n_controls = length(self.B);
        
        if nargin < 3 || isempty(st_labels)
            st_labels = char('0' + (1:dim).');
        end

        if nargin < 4 || isempty(c_labels)
            c_labels = char('0' + (1:n_controls).');
        end
        
        if length(st_labels) ~= dim
            error('Number of state labels given does not match the Hilbert space dimension.')
        end
        self.state_labels = st_labels;
        
        if length(c_labels) ~= n_controls
            error('Number of control labels given does not match the number of controls.')
        end
        self.control_labels = c_labels;
    end
  end
end
