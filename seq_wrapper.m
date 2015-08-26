classdef seq_wrapper < matlab.mixin.Copyable
% Control sequence wrapper handle class.
% Takes care of import, export and integration, calls RWA.
% Does not know or care about the underlying physics beyond the Hamiltonians.

% Ville Bergholm 2012-2014
    
  properties
    desc = ''   % system description
    TU          % time unit (s)
    dim         % dimension vector
    calibration % which calibration shall we use in import/export (one per carrier)
    coda = ''   % filename extension for using special calibrations

    tau        % pulse durations (TU)
    cum_tau    % cumsum(tau), final bin == inf
    amp, phi   % control fields, size(amp) == [n_bins, n_controls]. amp is in (1/TU)
    carrier    % carrier frequency for each control (1/TU)
    t0 = 0     % starting time for the sequence (for rotating terms etc.)
    
    H          % drift Hamiltonian (\hbar / TU)
    H0         % rotating frame generator Hamiltonian
    C          % control Hamiltonians

    %% flags
    flags = struct('use_LO', 0,...      % include a LO leak?
                   'safety', 500,...    % safety limit for eliminating rotating terms
                   'verbosity', 1,...   % should we print lots of stuff or not?
                   'in_florian', 0,...  % use Florian's field ordering during import
                   'in_reverse', 0)     % reverse the sequence during import
                   
    %% local oscillator leak properties
    omega_rabi_LO
    phi_LO

    %% rotating frame stuff
    H_rwa      % drift Hamiltonian (static terms)
    Hr_rwa     % drift Hamiltonian (rotating terms)
    C_rwa      % control H:s (static)
    Cr_rwa     % control H:s (rotating)
    amp_min    % minimum control amplitude at which each Cr_rwa operator remains significant
    amp_max    % maximum control amplitude that keeps fast modes negligible
    omega_rot  % rotation speeds
  
    amp_keep = [] % keep all Cr_rwa terms whose amp_min is below this
  end

  
  methods (Access = private)
  end
  
  methods (Static)
    function print_row(a_min, omega_rot)
    % Utility.
        fprintf('a_min:');
        fprintf(' %g,', a_min);
        fprintf('\nomega_rot:');
        fprintf(' %g,', omega_rot);
        fprintf('\n');
    end
  end
  
  methods
    function self = seq_wrapper(TU, H, H0, dim, C, carrier, flags, ord)
    % Constructor. Set up the Hamiltonians etc.
    % Carriers, if given, are normal frequencies (not angular).

        %% apply flags
        if nargin >= 7
            temp = fieldnames(flags);
            for k=1:length(temp)
                self.flags = setfield(self.flags, temp{k}, getfield(flags, temp{k}));
            end
        end

        %% reorder subsystems?
        if nargin >= 8
            % abusing state class
            H0 = state(H0, dim).reorder(ord).data;
            H  = state(H,  dim).reorder(ord).data;
            for k=1:length(C)
                C{k} = state(C{k},  dim).reorder(ord).data;
            end
            
            % and finally dim
            dim = dim(ord);
        end

        self.TU = TU;
        self.dim = dim;
        self.H = H;
        self.set_H0(H0);  % define the rotating frame
        

        %% carrier freqs/controls
        if ischar(carrier)
            % carrier is a string defining the carrier freqs
            switch carrier
              case 'peaks'
                % NOTE this is roughly how Florian determines the resonance freqs from the spectrum of H.
                % The difference from the bare H_e resonance freqs is
                % only about 3.5 kHz at theta_th, B0 = 20 Gauss, or
                % about 20 kHz at theta = 2.0, B0 = 30 Gauss. It gets
                % larger when theta = pi/2 though.
                fprintf('\nODMR peak carriers.')
                d = eig(self.H);
                temp = sum(d(1:2)) / 2;  % normally not resolved, so average them
                lo   = sum(d(3:4)) / 2;  % average of hyperfine peaks
                hi   = sum(d(5:6)) / 2;  % same
                carrier = [lo-temp, hi-temp] / (2*pi);
                
              case 'hyperfine'
                % For the rough SWAP gates.
                fprintf('\nHyperfine selective carriers.')
                d = eig(self.H);
                temp = sum(d(1:2)) / 2;  % normally not resolved, so average them
                %carrier = [d(3)-temp, d(4)-temp, d(5)-temp, d(6)-temp];
                carrier = [d(4)-temp, d(6)-temp] / (2*pi);
                
              otherwise
                error('Unknown carrier name.')
            end
        end
        self.set_carriers(carrier, C);
    end

    function set_H0(self, H0)
    % Defines the rotating frame.

        self.H0 = H0;
      
        %% drift Hamiltonian, RWA
        [self.H_rwa, self.Hr_rwa, a_min, ~, omega_rot] =...
            RWA(self.H0, self.H -self.H0, 0, self.flags.verbosity, true, self.flags.safety);
    
        if self.flags.verbosity >= 1
            seq_wrapper.print_row(a_min, omega_rot)
        end
    end

    function set_carriers(self, freqs, C)
    % Sets the carrier frequencies, does the RWA for the control Hamiltonians.
    % set_H0 must have been called before.

        if self.flags.use_LO
          %% leaking LO params (as the last "control")
          C{end+1} = C{1}; % just copy the first one... TODO
          freqs{end+1} = [2.71e9 * self.TU]; % omega_LO
                                                       
          self.omega_rabi_LO = 7/85 * 10e6 * 2*pi * self.TU;
          %self.omega_rabi_LO = 10/65 * 10e6 * 2*pi * self.TU;
          self.phi_LO = pi/4;
          fprintf('Using the LO leak.\n');
        end

        % number of control Hamiltonians, == length(C)
        n = length(freqs);
        
        %% control Hamiltonians, RWA
        self.carrier = [];
        self.C       = {};
        self.C_rwa   = {};
        self.Cr_rwa  = {};
        self.amp_min = {};
        self.amp_max = [];
        self.omega_rot = {};

        % loop over control ops
        ind = 1;
        for j=1:n
          if self.flags.verbosity >= 1
            fprintf('\nControl operator %d:\n', j);
          end

          % list of carrier freqs for this control operator
          temp = freqs{j};
          for k=1:length(temp)
            self.carrier(ind) = temp(k);
            self.C{ind} = C{j};
            [self.C_rwa{ind}, self.Cr_rwa{ind}, self.amp_min{ind}, self.amp_max(ind), self.omega_rot{ind}] =...
                RWA(self.H0, C{j}, 2*pi * temp(k), self.flags.verbosity, true, self.flags.safety);
        
            if self.flags.verbosity >= 1
                seq_wrapper.print_row(self.amp_min{ind}, self.omega_rot{ind})
            end
            ind = ind+1;
          end
        end
        fprintf('\n\n');

        % keep the amp_keep information if possible
        if length(self.amp_keep) ~= length(self.carrier)
            self.amp_keep = inf(size(self.carrier));  % keep all
        end
    end

    function set_desc(self, desc, calibration)
    % Set system description and calibration data.

        self.desc = desc;
        self.calibration = calibration;
    end
    
    
    function [H, HC, dim, desc, HC_labels, control_type, control_par] = const_rwa(self)
    % Extract and return the constant RWA terms for constructing a
    % standard bilinear control system for DYNAMO etc.

        % Total drift Hamiltonian in the RWA
        H = self.H_rwa;

        % control Hamiltonians in the RWA (ignoring crosstalk and fast modes)
        C = self.C_rwa;

        dim = self.dim;
        
        n_carr = length(self.carrier);
        if self.flags.use_LO
            n_carr = n_carr -1; % the LO leak is not controllable
        end

        desc = 'RWA. Carriers:';
        for k=1:n_carr
            range = 2*k-1:2*k;
            HC(range) = {C{k}(0), C{k}(pi/2)}; % x and y phases

            desc = sprintf('%s, f_%d: %.8g MHz', desc, k, self.carrier(k) / self.TU * 1e-6);
            HC_labels(range) = {sprintf('x%d', k), sprintf('y%d', k)};

            control_type(range) = 'mm';
            temp = self.amp_max(k);
            control_par(range) = {[-1, 2] * temp, [-1, 2] * temp};
        end

        % add system description
        desc = sprintf('%s\n%s', desc, self.desc);
    end
    

    function post_import(self)
    % Does some basic housekeeping after changing the sequence.
        self.cum_tau = cumsum(self.tau);
        self.cum_tau(end) = inf;  % final bin goes on forever
        %fprintf('T = %.8f\n', sum(self.tau));
    end
    
    
    function import_better(self, tau, fields, in_SI, in_xy)
    % Imports a control sequence.

        % control fields
        n_bins = size(fields, 1);
        n_controls = size(fields, 2) / 2;
        n_carriers = length(self.carrier);

        self.amp = zeros(n_bins, n_carriers);
        self.phi = zeros(n_bins, n_carriers);

        % convert
        temp = 1:2:2*n_controls;
        if in_xy
            % in cartesian coordinates
            x = fields(:, temp);
            y = fields(:, temp+1);
            amp = sqrt(x.^2 +y.^2);
            phi = atan2(y, x);
        else
            % already polar
            amp = fields(:, temp);
            phi = fields(:, temp+1);
        end

        if in_SI
            % imported data is in SI units (tau (s), f (Hz))
            tau = tau / self.TU;
            % to angular freq
            amp = amp * 2*pi * self.TU;
        end
        
        % set
        self.tau = tau;
        self.amp(:, 1:n_controls) = amp;
        self.phi(:, 1:n_controls) = phi;
        
        if self.flags.use_LO
            % add the LO leak as the last "control"
            self.amp(:, n_controls+1) = self.omega_rabi_LO;
            self.phi(:, n_controls+1) = self.phi_LO;
            temp = n_carriers -1;
        else
            temp = n_carriers;
        end
        if n_controls ~= temp
            error('Wrong number of controls.')
        end

        %self.AWG_phase_bug();
        self.post_import();
    end
    
    
    function import(self, a)
    % Imports a control sequence array.

        % split it up
        self.import_better(a(:, 1), a(:, 2:end), true, false);
    end

    
    function AWG_phase_bug(self)
    % Simulates a possible bug in AWG code, resetting phases for each bin.

        disp('Simulating AWG phase bug.')
        temp = length(self.carrier);
        if self.flags.use_LO
            temp = temp -1;
        end
        % starting times for the bins
        t_start = cumsum(self.tau) -self.tau;
        % subtract the thus far accumulated phase
        for k=1:temp
            self.phi(:,k) = self.phi(:,k) -2*pi * self.carrier(k) * t_start;
        end
        % normalize the angles
        temp = self.phi / (2 * pi);
        temp = temp -floor(temp); % remove full rotations
        self.phi = 2*pi * temp;
    end
    
    
    function [out, desc] = import_file(self, fname, just_array)
    % Imports an NV- control sequence from a file.

        if nargin < 3
            just_array = false;
        end
        
        desc = [];
        id = fopen(fname, 'r', 'ieee-be', 'UTF-8');
        % discard comment lines until we hit an empty one
        while 1
            s = fgetl(id);
            desc = sprintf('%s%s\n', desc, s);
            disp(s);
            if strcmp(s, '')
                break
            end
        end
        % read in an unknown number of columns (max 20), fill the extra cols with NaN
        out = textscan(id, [repmat('%f', 1, 20), '%*[^\n]'], 'CollectOutput',true);
        out = out{1};
        % fucking hell! find the last not-NaN entry
        out = out(:, 1:find(~isnan(out(1,:)), 1, 'last'));
        fclose(id);
        rows_cols = size(out)
        
        if self.flags.in_florian
            % florian's ordering: first phases, then rabi f:s
            ind = 2:2:size(out, 2);
            temp = out(:,ind);
            out(:,ind) = out(:,ind+1);
            out(:,ind+1) = temp;
            out(:,1) = out(:,1) / 1.2e9;  % samples to s
            disp('Florian''s order!!!')
        end
        if self.flags.in_reverse
            % time-reverse the sequence
            out = flipud(out);
            disp('Reversing the sequence!!!')
        end

        % FIXME heuristic: amp value or Rabi frequency?
        temp = max(max(out(:, 2:2:end)));
        if temp < 2
            disp('assuming amp values');
            % amp to rabi
            for k=1:length(self.carrier)
                ind = 2*k;
                out(:,ind) = rabi_to_amp(self.calibration{k}, out(:,ind), true, self.coda);
            end
        end
        
        if just_array
            return
        end
        self.import_better(out(:, 1), out(:, 2:end), true, false);
    end

    
    function out = export(self, filename, to_amp, desc)
    % Exports an NV- control sequence into a file.
    % The output array rows correspond to pulses/bins, each row
    % consisting of [tau (s), f_rabi_1 (Hz), \phi_1, f_rabi_2 (Hz), \phi_2, ...]
    %
    % Alternatively, instead of driving Rabi frequencies, the
    % driving amplitudes can be expressed as amplifier settings if
    % a calibration file is available.

        n_bins = size(self.amp, 1);
        n_controls = size(self.amp, 2);
        if self.flags.use_LO
            n_controls = n_controls-1;
        end
        
        % f_rabi to amp setting conversion
        if nargin < 4
            desc = 'no description\n';
            if nargin < 3
                to_amp = false;
            end
        end

        % from angular to plain frequency, in Hz
        temp = self.amp(:, 1:n_controls) / (2*pi) / self.TU;

        if to_amp
            % each carrier frequency has its calibration in a separate file
            for k=1:length(self.carrier)
                temp(:, k) = rabi_to_amp(self.calibration{k}, temp(:, k), false, self.coda);
            end
        end

        out = [self.tau * self.TU];
        % splice the amp and phi tables together
        for k=1:n_controls
            ind = 2*k;
            out(:, ind)   = temp(:, k);
            out(:, ind+1) = self.phi(:, k);
        end
        
        if nargin < 2
            % just return the array, don't create a file
            return
        end

        % write the file
        filename = [filename, '.txt']
        id = fopen(filename, 'w', 'ieee-be', 'UTF-8');
        fprintf(id, 'Dynamo control sequence\n');
        %fprintf(id, 'Generated %s\n', d.config.date_UTC);
        fprintf(id, '=======================\n');
        fprintf(id, '%s\n', desc);
        %fprintf(id, 'Error: %g\n', d.config.error_func(d));
        fprintf(id, 'Rows correspond to pulses, each row consists of\n');
        if to_amp
            fprintf(id, '[tau (s), amp_1, phi_1, amp_2, phi_2, ...]\n\n');
        else
            fprintf(id, '[tau (s), f_rabi_1 (Hz), phi_1, f_rabi_2 (Hz), phi_2, ...]\n\n');
        end
        for k=1:n_bins
            fprintf(id, '% 12.6g ', out(k,:));
            fprintf(id, '\n');
        end
        fclose(id);
    end

    
    function [u, t] = get_control(self, t, c)
    % Control signal corresponding to the c:th carrier freq.
    % Sampled at times in t, u(t) = f(t)*exp(1i*phi(t))

        if nargin < 3
            % all controls
            c = 1:length(self.carrier);
            
            if nargin < 2 || isempty(t)
                % whole sequence, linearly sampled
                t = linspace(0, sum(self.tau), 2*length(self.tau));
            end
        end

        u = zeros(length(t), length(c));
        for k=1:length(t)
            tt = t(k);
            % find the correct bin
            bin = find(tt <= self.cum_tau, 1);
            u(k, :) = self.amp(bin, c) .* exp(1i * self.phi(bin, c));
            %u(k, :) = self.amp(bin, c) .* cos(self.phi(bin, c));
        end
    end

    

    function H = get_H(self, t)
    % Full time-dependent Hamiltonian
        
        % find the correct bin
        bin = find(t <= self.cum_tau, 1);
    
        t = t +self.t0;
        
        % drift
        H = self.H;
        % loop over control fields
        for k=1:length(self.carrier)
            H = H +self.amp(bin,k) * cos(2*pi * self.carrier(k) * t +self.phi(bin,k)) * self.C{k};
        end
    end

    
    function H = get_H_rwa(self, t)
    % RWA Hamiltonian
    
        % find the correct bin
        bin = find(t <= self.cum_tau, 1);
  
        t = t +self.t0;
        
        %% drift
        % static terms
        H = self.H_rwa;
        
        % the most significant rotating terms (no phi)
        for j=1:length(self.Hr_rwa)
            H = H +self.Hr_rwa{j}(t);
        end
        
        %% loop over control fields
        for k=1:length(self.carrier)
            phi = self.phi(bin,k);
            
            % static component
            temp = self.C_rwa{k}(phi);
            
            % the most significant rotating components
            for j=1:length(self.Cr_rwa{k})
                temp = temp + self.Cr_rwa{k}{j}(t, phi);
            end
            H = H +self.amp(bin,k) * temp;
        end
    end

    
    function U = int(self, t, in_labframe)
    % Integrate the sequence to the time instance(s) in t, return
    % the propagator(s) in the lab frame.
    % This "exact" integrator tracks every oscillation of the
    % control fields and is extremely slow. It is meant to be used
    % as a benchmark only.
 
        if nargin < 3
            in_labframe = false;
        end
        
        U = self.do_int(t, @(x) get_H(self, x));
    
        % convert to rotating frame?
        if ~in_labframe
            U = self.to_lab(t, U, true);
        end
    end


    function U = int_bin_rwa(self, bin)
    % Integrates a propagator for a single bin, in rotating frame.
        
        t_end   = sum(self.tau(1:bin));
        t_start = t_end -self.tau(bin);
        
        U = int_rwa(self, t_end, false, t_start);
        % pick only the one we care about
        U = U{end};
    end


    function U = int_rwa(self, t, in_labframe, t_start)
    % Integrate the sequence using RWA, results are in rotframe by default.
    % Integration is from t_start (or 0, if it is not given) to the points in the vector t.
    % Returns a cell vector of propagators.
        
        function rrr_ret = H_rot(rrr_t)
        % Hamiltonian in the rotating frame.
        % The rrr_ prefix is there to not fuck up the variables in the containing scope!

            rrr_ret = rrr_H; % start with the static part
            rrr_t = rrr_t +self.t0;
            
            % rotating terms of drift
            for rrr_j = 1:length(self.Hr_rwa)
                rrr_ret = rrr_ret +self.Hr_rwa{rrr_j}(rrr_t);
            end
            
            % rotating control terms
            for rrr_k = 1:length(self.carrier)
                temp = 0;
                for rrr_j = 1:rrr_keep(rrr_k)
                    temp = temp + rrr_Cr{rrr_k}{rrr_j}(rrr_t, rrr_phi(rrr_k));
                end
                rrr_ret = rrr_ret +rrr_amp(rrr_k) * temp;
            end
        end

        % parent function starts here
        if nargin < 4
            t_start = 0;
            if nargin < 3
                in_labframe = false;
                if nargin < 2
                    t = [];
                end
            end
        end
        if isempty(t)
            % integrate the whole sequence
            t = sum(self.tau);
        end
        if  t_start < 0
            error('t_start cannot be negative.')
        end
        if any(t < t_start)
            error('t values must be larger than (or equal to) t_start.')
        end
        
        % floating point tolerance (err on the side of faster computation)
        tol = 1e-10;
        
        % make a copy
        t_todo = t;

        if 0
            U = self.do_int(t, @(x) get_H_rwa(self, x));
        else
        % TEST: bin by bin, faster. 
        % This is bullshit. It seems that self.x references are slow as hell.
        rrr_Cr = self.Cr_rwa;  % shorthand for use in H_rot
        
        done = false;
        U_last = eye(size(self.H0));
        U = {};
        
        % loop over the bins
        for bin = 1:length(self.tau)
            % end of the bin (inf for the last bin)
            t_bin_end = self.cum_tau(bin);

            if t_bin_end <= t_start +tol
                % we aren't interested in this bin
                continue
            end
            
            % find the t:s belonging to this bin
            temp = (t_todo <= t_bin_end +tol);
            t_in_this_bin = t_todo(temp);
            % time points in the remaining bins
            t_todo = t_todo(~temp);

            % required time points for the integrator
            if isempty(t_todo)
                done = true;
                t_int = [t_start];
            else
                % for the next bin we need the full propagator for
                % this bin, hence t_bin_end is required
                t_int = [t_start, t_bin_end];
            end
            % add the points we are interested in, sort
            t_int = unique(union(t_in_this_bin, t_int));
            
            % fucking ode45 needs at least three distinct points
            % (two are guaranteed by the above code (unless initially t = [0]))
            if length(t_int) < 3
                % insert a superfluous inside point
                temp = t_int(1);
                t_int = [temp, 0.5*(temp +t_int(2)), t_int(2)];
            end

            rrr_phi = self.phi(bin,:);
            rrr_amp = self.amp(bin,:);

            % build the static part of the bin Hamiltonian
            rrr_H = self.H_rwa;
            % loop over controls
            for k = 1:length(self.carrier)
                % static component
                rrr_H = rrr_H +rrr_amp(k) * self.C_rwa{k}(rrr_phi(k));
            
                % number of rotating terms to keep for each
                % carrier, depending on the amp value in this bin
                %rrr_keep(k) = length(find(amp(k) >= self.amp_min{k}));
                %rrr_keep(k) = length(self.amp_min{k}); % keep'em all
                rrr_keep(k) = length(find(self.amp_min{k} <= self.amp_keep(k)));
            end

            % integrate the propagators
            temp = self.do_int(t_int, @H_rot, U_last);
            % pick only the ones we care about
            U = horzcat(U, temp(ismember(t_int, t_in_this_bin)));
            
            if done
                break
            end
            % this bin is done
            U_last = temp{end};
            t_start = t_bin_end;
        end
        end
        
        % convert to lab frame
        if in_labframe
            U = self.to_lab(t, U);
        end
    end
    
    
    function U = int_expm(self, t, in_labframe)
    % Integrate the sequence using RWA, keeping only static
    % terms. This means that expm is enough.
    % Results are in rotframe by default.
        
        if nargin < 3
            in_labframe = false;
            
            if nargin < 2
                % integrate the whole sequence
                t = [];
            end
        end
        if isempty(t)
            % integrate the whole sequence
            t = sum(self.tau);
        end

        
        U_last = eye(size(self.H0));
        t_start = 0;
        U = {};
        
        % loop over bins
        for bin = 1:length(self.tau) 
            % end of the bin
            t_bin_end = self.cum_tau(bin);
            
            % find the t:s belonging to this bin
            t_in_this_bin = t(t_start <= t & t < t_bin_end);
            
            % bin Hamiltonian
            H = self.H_rwa;
            % loop over control fields
            for k=1:length(self.carrier)
                % static component
                H = H +self.amp(bin,k) * self.C_rwa{k}(self.phi(bin,k));
            end

            % propagators
            for s = t_in_this_bin
                U = horzcat(U, expm(-1i * (s-t_start) * H) * U_last);
            end
            
            % this bin is done
            U_last = expm(-1i * self.tau(bin) * H) * U_last;
            t_start = t_bin_end;
        end
        
        % convert to lab frame
        if in_labframe
            U = self.to_lab(t, U);
        end
    end

    
    
    function [ret, t_out] = do_int(self, t, H, U_initial)
    % integrate the system propagator, t is a vector of increasing time instances

        odeopts = odeset('RelTol', 1e-8, 'AbsTol', 1e-8); %, 'MaxStep', aaa);

        if nargin < 4
            temp = H(0);
            U_initial = eye(size(temp));
        end
        y0 = vec(U_initial);

        function ret_ = odefun(t_, y_)
            ret_ = vec(-1i * H(t_) * inv_vec(y_));
        end
        
        [t_out, y_out] = ode45(@odefun, t, y0, odeopts);
        for k = 1:length(t_out)
            ret{k} = inv_vec(y_out(k, :));
        end
    end

  
    function U = to_lab(self, t, U, invert)
    % rotating frame propagator U(t <- 0) to lab frame.
    % If invert is true, will do the opposite transformation.

        if nargin < 4
            invert = false;
        end

        % take sequence starting time into account
        t0 = self.t0;
        t = t +t0;

        % invert the transformation by negating the times
        if invert
            t = -t;
            t0 = -t0;
        end

        V = expm(1i * self.H0 * t0);
        for k=1:length(t)
            U{k} = expm(-1i * self.H0 * t(k)) * U{k} * V;
        end
    end
  end
end
