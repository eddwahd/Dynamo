classdef dynamo < matlab.mixin.Copyable
% Copyable handle class for DYNAMO optimizer objects.
%
% Contains the optimization task, system description, control
% sequence, various options and statistics. This class glues together
% the functionality of the cache, qsystem and control classes.
%
% Governing equation: \dot(X)(t) = -(A +\sum_k u_k(t) B_k) X(t) = -G(t) X(t)
    
% Shai Machnes   2010-2011
% Ville Bergholm 2011-2012


  properties
    config   % configuration information, other metadata
    system   % description of the physical system
    seq      % control sequence
    opt      % optimization options
    stats    % optimization statistics
  end

  properties (Transient)
    cache % Do not save the cache on disk since it may be huge and can always be recomputed.
  end

  methods (Static)
    function ret = version()
    % Returns the current DYNAMO version.
        ret = '1.3 alpha14';
    end


    function obj = loadobj(obj)
    % Re-initializes the cache (which is not saved) during loading.
        if isstruct(obj)
            error('Backwards compatibility of saved objects not yet implemented.')
        end
        obj.cache_init();
    end
  end

  methods (Access = protected)
    function cp = copyElement(self)
    % Override the default copyElement method to provide deep copies.
    % This is necessary since some of the data members are Copyable handle classes themselves.
        
        % Make a shallow copy of everything
        cp = copyElement@matlab.mixin.Copyable(self);
        % Make a deep copy of all handle-type properties
        cp.system = copy(self.system);
        cp.seq    = copy(self.seq);
        cp.cache  = copy(self.cache);
    end
  end
  
  methods
    function self = dynamo(task, initial, final, H_drift, H_ctrl, L_drift)
    % Constructor

        if nargin < 6
            L_drift = 0;
        end
        
        task = lower(task);

        %% Some basic data provenance

        config.version = dynamo.version();
        % Local time. TODO UTC or local time with timezone specifier would be better, but apparently MATLAB doesn't do that.
        config.date = datestr(now(), 31);
        config.task = task;
        config.UL_hack = false;
        
        [system_str, rem] = strtok(task);
        [task_str, rem] = strtok(rem);
        [extra_str, rem] = strtok(rem);
          

        %% Description of the physical system

        sys = qsystem();

        % TODO FIXME temporary fix: sparse to full
        L_drift = full(L_drift);
        H_drift = full(H_drift);
        for k = 1:length(H_ctrl)
            H_ctrl{k} = full(H_ctrl{k});
        end
        input_dim = [size(initial, 2), size(final, 2)]; % check the validity of the inputs

        out = 'Target operation:';
        switch system_str
          case {'s'}
            %% Closed system S
            if nargin == 6
                error('L_drift not used in closed systems.')
            end
            % the generator is always Hermitian and thus normal => use exact gradient

            sys.hilbert_representation(initial, final, H_drift, H_ctrl);

            switch task_str
              case 'state'
                % TEST more efficient Hilbert space implementation
                out = strcat(out, ' mixed state transfer');
                if any(input_dim == 1)
                    error('Initial and final states should be state operators.')
                end
                config.error_func = @error_real;
                config.f_max = 0.5 * (1 + norm2(sys.X_initial) / norm2(sys.X_final));
                config.gradient_func = @gradient_g_mixed_exact;
                config.UL_hack = true;
                
              case {'ket', 'gate'}
                if strcmp(task_str, 'ket')
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
                
                if strcmp(extra_str, 'phase')
                    out = strcat(out, ' (with global phase (NOTE: unphysical!))');
                    config.error_func = @error_real;
                else
                    out = strcat(out, ' (ignoring global phase)');
                    config.error_func = @error_abs;
                end
                config.f_max = 1;
                config.gradient_func = @gradient_g_exact;
                
              otherwise
                error('Unknown task.')
            end
            out = strcat(out, ' in a closed system.\n');


          case {'sb'}
            %% Open system S with bath B
            % The generator isn't usually normal, so we cannot use the exact gradient method

            switch task_str
              case 'state'
                out = strcat(out, ' quantum state transfer');
                sys.vec_representation(initial, final, H_drift, L_drift, H_ctrl);

              case 'gate'
                out = strcat(out, ' quantum gate');
                if any(input_dim == 1)
                    error('Initial and final states should be unitary operators.')
                end
                sys.vec_gate_representation(initial, final, H_drift, L_drift, H_ctrl);
        
              otherwise
                % TODO arbitrary quantum maps
                error('Unknown task.')
            end

            if strcmp(extra_str, 'overlap')
                % overlap error function
                % NOTE simpler error function and gradient, but final state needs to be pure
                % TODO should not be used with gate task
                out = strcat(out, ' (overlap)');
                config.error_func = @error_real;
                config.gradient_func = @gradient_g_1st_order;
                config.f_max = 1;
            else
                % full distance error function
                config.error_func = @error_open;
                config.gradient_func = @gradient_open_1st_order;
            end
            out = strcat(out, ' in an open system under Markovian noise.\n');


          case {'se'}
            %% Closed system S + environment E
            if nargin == 6
                error('L_drift not used in closed systems.')
            end

            switch task_str
              case {'ket', 'state'}
                error('Not implemented yet.')

              case 'gate'
                out = strcat(out, ' unitary gate on S');
                if any(input_dim == 1)
                    error('Initial and final states should be unitary operators.')
                end

                M = input_dim(1) / input_dim(2);
                if floor(M) ~= M
                    error('Initial state must be a unitary on S+E, final state a unitary on S.');
                end
                sys.hilbert_representation(initial, kron(final, eye(M)), H_drift, H_ctrl);

                config.dimS = input_dim(2);
                config.error_func = @error_partial;
                config.gradient_func = @gradient_partial_exact;
                
              otherwise
                error('Unknown task.')
            end
            out = strcat(out, ' in a closed system S+E.\n');            


          case {'seb'}
            %% Open system S + environment E with bath B
            error('Not implemented yet.')
    
          otherwise
            error('Unknown system specification.')
        end
        fprintf(out);
        if sys.liouville
            fprintf('Liouville');
        else
            fprintf('Hilbert');
        end
        fprintf(' space dimension: %d\n\n', length(sys.X_final));
          
        % Calculate the squared norm |X_final|^2 to scale the fidelities with.
        % We use the Hilbert-Schmidt inner product (and the induced Frobenius norm) throughout the code.
        sys.norm2 = norm2(sys.X_final);
        
        % store the prepared fields
        self.config = config;
        self.system = sys;

        % init miscellaneous things
        self.opt.ui_fig = [];
        self.opt.max_violation = 0; % track the worst gradient approximation violation
    end


    function cache_init(self)
    % Set up cache after the number of time slots changes.
    % This is where all the bad code went.
        
        % error_open needs a full reverse propagator.
        if isequal(self.config.error_func, @error_open)
            L_end = eye(length(self.system.X_final)); % L: full reverse propagator
        else
            L_end = self.system.X_final'; % L: X_final' propagated backwards
        end

        % exact gradient? we need the eigendecomposition data.
        use_eig = false;
        temp = self.config.gradient_func;
        if isequal(temp, @gradient_g_exact)...
                || isequal(temp, @gradient_g_mixed_exact)...
                || isequal(temp, @gradient_partial_exact)
            use_eig = true;
        end

        % UL_hack: mixed states in a closed system
        self.cache = cache(self.seq.n_timeslots(), self.system.X_initial, L_end, use_eig, self.config.UL_hack);
    end


    function seq_init(self, n_timeslots, tau_par, varargin)
    % Create the control sequence and a matching cache.
    % The varargin are the control_type and control_par cell vectors.

        n_controls = length(self.system.B);
        self.seq = control_seq(n_timeslots, n_controls, tau_par, varargin{:});
        if any(self.seq.control_type(~self.system.B_is_Hamiltonian) == '.')
            disp('Warning: Liouvillian control ops with possibly negative control values.')
        end
    
        self.cache_init();
    end


    function mask = full_mask(self, optimize_tau)
    % Returns a full control mask.
        
        if nargin < 2
            optimize_tau = false;
        end
        n_timeslots = self.seq.n_timeslots();
        n_controls = self.seq.n_controls();

        % shape vectors
        f_shape = [n_timeslots, n_controls];
        t_shape = [n_timeslots, 1];

        %% Build the control mask

        fprintf('Tau values ');
        if optimize_tau
            fprintf('optimized.\n');
            mask = [true(f_shape), true(t_shape)];
        else
            fprintf('fixed.\n')
            mask = [true(f_shape), false(t_shape)];
        end
    end
    
    
    function update_controls(self, x, control_mask)
    % Updates selected controls.
    %
    %  x: vector of control values, control_mask: corresponding mask.
    %
    %  Updates control values for which control_mask is true.
    %  Makes the changed timeslots and stuff that depends on them stale.

        old = self.seq.get();
        
        if nargin < 3
            control_mask = true(size(old)); % full mask
        end

        % make a trial copy of the new controls
        new = old;
        new(control_mask) = x;

        % see which timeslots have changed
        changed_t_mask = any(new ~= old, 2);

        if any(changed_t_mask)
            % actually update the controls
            self.seq.set(new);
            self.cache.mark_as_stale(changed_t_mask);
        end
    end

    
    function cache_refresh(self)
    % Performs all the queued computations using the cache subsystem.
        self.cache.refresh(self.system, self.seq.tau, self.seq.fields);
    end
    

    function cache_fill(self)
    % Will invalidate everything, then re-calc everything in the cache.
    % Used mostly for debugging (since it essentially overrides all matrix-op optimization mechanisms).

        self.cache.invalidate();

        self.cache.H_needed_now(:) = true;
        self.cache.P_needed_now(:) = true;
        self.cache.U_needed_now(:) = true;
        self.cache.L_needed_now(:) = true;
        self.cache.g_needed_now    = true;
        
        self.cache_refresh();
    end


    function ret = g_func(self)
    % Computes the auxiliary function g := trace(X_f^\dagger * X(t_n)).
    % Used both for the goal function as well as its gradient.
        
        self.cache.g_needed_now = true;
        self.cache_refresh();
        ret = self.cache.g;
    end
  
  
    function ret = X(self, k)
    % Returns X(t_k), the controlled system at time t_k.
    % If no k is given, returns the final state X(t_n).

        if nargin < 2
            k = length(self.cache.H);
        end
        % U{k+1} is the system at t = sum(tau(1:k)) = t_{k}
        self.cache.U_needed_now(k+1) = true;
        self.cache_refresh();
        ret = self.cache.U{k+1};
    end


    function plot_seq(self, varargin)
    % Plots the control sequence. Wrapper.

        ax = self.seq.plot(self.system.control_labels, varargin{:});
        title(ax, self.system.description);
    end


    function plot_stats(self, ax0)
    % Plots optimization stats.
        
        [ax, h1, h2] = plotyy(ax0, self.stats.wall_time, abs(self.stats.error), ...
                              self.stats.wall_time, self.stats.integral, 'semilogy', 'plot');
        % For some strange reason only the first set of axes inherits the
        % parent axes' properties, so we need to set the plotstyle separately for both.
        set_plotstyle(ax(1));
        set_plotstyle(ax(2));
        title('Optimization Statistics')
        xlabel('wall time (s)')
        ylabel(ax(1), 'normalized error')
        ylabel(ax(2), 'control integral')
        grid on
        set(h2, 'LineStyle','--')
    end


    function plot_X(self, ax, full, dt)
    % Plots the evolution of the initial system as a function of time.
    % TODO for now it only handles kets and state ops

        n = self.seq.n_timeslots();
        
        if nargin < 2
            ax = gca();
        end
        if nargin < 3
            full = true;
        end

        if full
            % things that don't change and aren't deleted by cla
            set_plotstyle(ax);
            title(ax, self.system.description);
            xlabel(ax, 'time');
            ylabel(ax, 'probability');
            grid(ax, 'on')
            set(ax, 'NextPlot','replacechildren'); % so plot() won't reset these properties
        else
            %cla(ax);
        end

        % what should we plot?
        if self.system.liouville
            state_probs = @prob_stateop;
        else
            state_probs = @prob_ket;
        end
        
        if nargin < 4
            % one plot point per timeslot
            for k = 0:n
                res(k+1, :) = state_probs(self.X(k));
            end
            t = [0; cumsum(self.seq.tau)];
        else
            % use the given dt for plotting
            t = [0];
            X = self.X(0);
            res(1, :) = state_probs(X);

            for k = 1:n
                X_end = self.X(k); % make sure the cache is up-to-date
                G = self.cache.H{k};
                tau = self.seq.tau(k);
                
                n_steps = ceil(tau / dt); % at least one point per timeslot
                step = tau / n_steps;
                P = expm(-step * G);
                for q = 1:n_steps
                    X = P * X;
                    res(end+1, :) = state_probs(X);
                end
                temp = t(end);
                t = [t, linspace(temp+step, temp+tau, n_steps)];
                X = X_end; % stability...
            end
        end
        plot(ax, t, res);
        axis(ax, [0, t(end), 0, 1]);
        legend(self.system.state_labels);


        function p = prob_stateop(x)
        % Returns the diagonal of a vectorized state operator.
        % NOTE: due to the horrible scoping rules of MATLAB, we use small x
        % here as not to nuke the capital X in the parent function.
            p = real(diag(inv_vec(x)));
        end
    
        function p = prob_ket(x)
        % Returns the absolute values squared of ket vector elements.
        % NOTE: due to the horrible scoping rules of MATLAB, we use small x
        % here as not to nuke the capital X in the parent function.
            p = real(x .* conj(x));
        end
    end
  end
end
