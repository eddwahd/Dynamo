classdef qsystem
% Defines a quantum system.

% Ville Bergholm 2011-2012
    
  properties
    liouvillian      % Are the objects below in Liouville space or Hilbert space?
    labels = {}      % Hilbert space state labels
    A                % drift generator
    B                % cell vector of control generators
    B_is_Hamiltonian % true if corresponding B item is a Hamiltonian
    M
    X_initial        % initial state
    X_final          % final state
    norm2            % squared norm of final state
  end
    
  methods
    function self = std_representation(self, i, f)
    % X_* are Hilbert space objects (kets or operators)
        self.X_initial = i;
        self.X_final   = f;
    end


    function self = vec_representation(self, i, f)
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


    function self = vec_gate_representation(self, i, f)
    % X_* are Liouville space operators corresponding to
    % vec-torized unitary gates.
        self.X_initial = lrmul(i, i'); % == kron(conj(i), i);
        self.X_final   = lrmul(f, f'); % == kron(conj(f), f);
    end


    function self = hilbert(self, H_drift, H_ctrl)
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


    function self = liouville(self, H_drift, L_drift, H_ctrl)
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


    function self = set_labels(self, labels)
    % Label the states. labels is a cell vector of strings.

        if length(labels) ~= self.dimension()
            error('Number of labels given does not match the Hilbert space dimension.')
        end
        self.labels = labels;
    end
  end
end
