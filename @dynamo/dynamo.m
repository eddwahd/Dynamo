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
    config
    system
    seq
    opt
    stats
  end

  properties (Transient)
    cache
  end

  methods (Static)
    function ret = version()
    % Returns the current DYNAMO version.
        ret = '1.3 alpha11';
    end
  end
  
  methods
    function self = dynamo(task, initial, final, H_drift, H_ctrl, L_drift)
    % Constructor

        task = lower(task);

        %% Some basic data provenance

        config.version = dynamo.version();
        % Local time. TODO UTC or local time with timezone specifier would be better, but apparently MATLAB doesn't do that.
        config.date = datestr(now(), 31);
        config.task = task;
        config.L_is_propagator = false;
        
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
            if nargin ~= 5
                error('Too many parameters.')
            end
    
            switch task_str
              case 'state'
                out = strcat(out, ' mixed state transfer');
                % TODO more efficient Hilbert space implementation?
                sys = sys.vec_representation(initial, final);
                sys = sys.liouville(H_drift, 0, H_ctrl);
                % g is always real, positive in this case so error_abs would work just as well
                config.error_func = @error_real;
        
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

                sys = sys.std_representation(initial, final);
                sys = sys.hilbert(H_drift, H_ctrl);
        
                if strcmp(extra_str, 'phase')
                    out = strcat(out, ' (with global phase (NOTE: unphysical!))');
                    config.error_func = @error_real;
                else
                    out = strcat(out, ' (ignoring global phase)');
                    config.error_func = @error_abs;
                end
        
              otherwise
                error('Unknown task.')
            end
            out = strcat(out, ' in a closed system.\n');
            
            % global maximum of the quality function f_max
            config.f_max = sqrt(norm2(sys.X_initial) / norm2(sys.X_final));

            % the generator is always Hermitian and thus normal => use exact gradient
            config.gradient_func = @gradient_exact;

    
          case {'sb'}
            %% Open system S with bath B
            switch task_str
              case 'state'
                out = strcat(out, ' quantum state transfer');
                sys = sys.vec_representation(initial, final);

              case 'gate'
                out = strcat(out, ' quantum map');
                if any(input_dim == 1)
                    error('Initial and final states should be operators.')
                end
                sys = sys.vec_gate_representation(initial, final);
        
              otherwise
                error('Unknown task.')
            end
            sys = sys.liouville(H_drift, L_drift, H_ctrl);

            % The generator isn't usually normal, so we cannot use the exact gradient method
            self.opt.max_violation = 0; % track the worst violation

            if strcmp(extra_str, 'overlap')
                % TEST, simple overlap goal function
                out = strcat(out, ' (overlap)');
                config.error_func = @error_real;
                config.gradient_func = @gradient_first_order_aprox;
                config.f_max = 1;
            else
                config.error_func = @error_open;
                config.L_is_propagator = true; % L: full reverse propagator
            end
            out = strcat(out, ' in an open system under Markovian noise.\n');


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
        if sys.liouvillian
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
    end


    function cache_init(self, n_timeslots)
    % Set up cache after the number of time slots changes.
        
    % This is where all the bad code went.
        
        if self.config.L_is_propagator
            L_end = eye(length(self.system.X_final)); % L: full reverse propagator
        else
            L_end = self.system.X_final'; % L: X_final' propagated backwards
        end

        self.cache = cache(n_timeslots, self.system.X_initial, L_end, isequal(self.config.gradient_func, @gradient_exact));
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
            self.seq = self.seq.set(new);
            
            self.cache.H_is_stale(changed_t_mask) = true;
            self.cache.P_is_stale(changed_t_mask) = true;
    
            % Propagate the H_is_stale to the U and Ls.
            self.cache.U_is_stale( (find(self.cache.H_is_stale, 1, 'first')+1):end) = true;
            self.cache.L_is_stale(1:find(self.cache.H_is_stale, 1, 'last'))         = true;
    
            self.cache.g_is_stale = true;
            self.cache.g = NaN;
        end
    end

    
    function cache_refresh(self)
    % Performs all the queued computations using the cache subsystem.
        self.cache.refresh(self.system, self.seq.tau, self.seq.fields);
    end
    
    
    function ret = X(self, k)
    % Returns X(t_k), the controlled system at time t_k.
    % If no k is given, returns the final state X(t_n).
        
        if nargin < 2
            k = length(self.cache.H);
        end
        
        % U{k} is the system at t = sum(tau(1:(k-1))) = t_{k-1}
        self.cache.U_needed_now(k+1) = true;
        self.cache_refresh();
        ret = self.cache.U{k+1};
    end
  
    function ret = g_func(self, use_trace)
    % Computes the auxiliary function g := trace(X_f^\dagger * X(t_n)).
    % Used both for the goal function as well as its gradient.

        if nargin < 2
            use_trace = true; % default: apply trace
        end
        
        if use_trace && ~self.cache.g_is_stale 
            ret = self.cache.g;
            return;
        end

        % g can be computed using any slice k \in [1, n+1]: g = trace(L_k * U_k).
        % Try to figure out which k requires least additional computation.
        k = self.cache.g_setup_recalc();
        self.cache_refresh();

        if use_trace
            ret = trace_matmul(self.cache.L{k}, self.cache.U{k});
            self.cache.g = ret;
        else
            ret = self.cache.L{k} * self.cache.U{k};
            self.cache.g = trace(ret);
        end
        self.cache.g_is_stale = false;
    end
  end
end
