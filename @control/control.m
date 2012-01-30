classdef control
% Class for control sequences.

% Ville Bergholm 2011-2012

  properties
    tau_par
    control_type
    control_par

    raw
    tau
    tau_deriv
    fields
    fields_deriv
  end

  methods
    function self = control(n_timeslots, n_controls, tau_par, control_type, control_par)
    %  tau_par: [T_min, T_delta], either a single entry with the totals or one such entry for each time slot.

        %% Check control types
        if nargin < 4
            control_type = char(zeros(1, n_controls) + '.');
        elseif length(control_type) ~= n_controls
            error('Length of the control_type vector does not match the number of controls.')
        end

        fprintf('Timeslots: %d\nControls per slot: %d + tau\n\n', n_timeslots, n_controls);
        
        % TODO
        %if any(control_type(self.system.B_is_superop) == '.')
        %    disp('Warning: Liouvillian control ops with possibly negative control values.')
        %end

        %% Check control parameters
        % For now control_par doesn't have separate entries for each timeslot
        if nargin < 5
            control_par = {};
        end

        %% Check tau parameters
        % Is tau_par a total time or a vector of delta_t:s?
        len_params = size(tau_par, 1);
        if len_params == 1
            % divide the parameters uniformly into timeslots
            tau_par = ones(n_timeslots, 1) * (tau_par / n_timeslots);
        elseif len_params ~= n_timeslots
            error('Length of tau_par does not match the number of control timeslots given.')
        end

        self.tau_par = tau_par;
        self.control_type = control_type;
        self.control_par = control_par;
    end

    
    function ret = get(self, control_mask)
    % Returns the raw controls corresponding to the mask given, or all
    % of them if no mask is given.

        if nargin == 1
            ret = self.raw;
        else
            ret = self.raw(control_mask);
        end
    end

    
    function self = set(self, raw)
    % Transform and set the controls using a diagonal transformation function.
    %
    % raw: raw, untransformed control values, size(raw) == [n_timeslots, n_controls + 1].
    % Given raw controls r_k, computes the transformed controls c_k(r_k) and their derivatives.
    % Uses the (fixed) control parameters.

        %% Check the number of control fields. The last column are the tau controls.
        n_controls = length(self.control_type);
        if size(raw, 2) ~= n_controls + 1
            error('Number of controls given does not match the number of control fields.')
        end

        self.raw = raw;
        
        % Tau control is the last one.
        t_raw = raw(:, end);

        % Returned tau values should always be positive, and not too large
        % if we're using the 1st order gradient approximation.
        % Stretchy bins with min and max duration:  t_raw = 0 <=> min duration
        self.tau       = self.tau_par(:,1) +0.5 * self.tau_par(:,2) .* (1-cos(t_raw));
        self.tau_deriv = 0.5 * self.tau_par(:,2) .* sin(t_raw);

        temp = size(raw) - [0, 1];
        self.fields = zeros(temp);
        self.fields_deriv = zeros(temp);

        for k=1:n_controls
            switch self.control_type(k)
              case '.'  % no transformation
                self.fields(:, k) = raw(:, k);
                self.fields_deriv(:, k) = 1;
        
              case 'p'  % strictly nonnegative, u_k = r_k^2
                self.fields(:, k) = raw(:, k).^2;
                self.fields_deriv(:, k) = 2 * raw(:, k);
                
              case 'm'  % minimum and delta, u_k = min + delta * 0.5 * (1 - cos(r_k))
                par = self.control_par{k};
                self.fields(:, k) = par(1) +par(2) * 0.5 * (1 - cos(raw(:, k)));
                self.fields_deriv(:, k) = par(2) * 0.5 * sin(raw(:, k));
      
              otherwise
                error('Unknown control type.')
            end
        end

        % TODO u_x = A*cos(phi), u_y = A*sin(phi) (not diagonal, but within the same bin, i.e. block diagonal)
        % fields_deriv(slot, c_j, raw_k) = d c_{sj} / d raw_{sk}
    end



    function ret = fluence(self, M)
    % Computes the total fluence of the control fields.
    %
    % The fluence is defined as F^2 = \int_0^T |H_c(t)|^2 dt,    H_c(t) = \sum_r H_r f_r(t).
    % If the control Hamiltonians are orthonormal, M_ij = (H_i, H_j) = \delta_{ij},
    % this gives F^2 = \sum_k \int_0^T |f_k(t)|^2 dt.
    %
    % For piecewise constant controls this reduces to F^2 = \sum_{ri} |a_{ri}|^2 \Delta t_{i}
    % F^2 has units of energy^2 * time.
    % TODO relation to signal energy? E = F^2 / \hbar?

        ret = 0;
        n_timeslots = size(self.fields, 1);
        for k = 1:n_timeslots
            c = self.fields(k, :); % all controls for this timeslot
            ret = ret +c * M * c.' * self.tau(k);
        end
    end

    
    
    function ret = integral(self)
    % Computes the time integral of the control fields.

        ret = transpose(self.tau) * self.fields;
    end
    


    function plot(self)
    % Plots a control sequence using superposed bars.

        % start times for pulses
        t = [0; cumsum(self.tau)];
        c = self.fields;
        nc = size(c, 2); % number of controls

        colormap(jet)
        set(gca, 'CLim', [0 nc])
        hold on
        for j=1:length(t)-1
            x = [t(j), t(j+1), t(j+1), t(j)];
    
            % use painter's algorithm
            [dummy, I] = sort(abs(c(j, :)), 2, 'descend'); % dummy instead of ~ for Octave
            for k=1:nc
                y = c(j, I(k));
                p(k) = patch(x, [0, 0, y, y], 1);
                set(p(k), 'FaceColor','flat',  'CData',I(k), 'CDataMapping','scaled')
            end
            if j == 1
                p_colors(I) = p; % HACK, to get the legend right
            end
        end
        axis tight
        title('Control sequence')
        xlabel('Time')
        ylabel('Control amplitude')
        grid on
        legend(p_colors, char('0' + (1:nc)'));
        hold off
        %stairs(t, c)
        %c = [c; zeros(1,nc)]; % final time step is a dummy
        %for j=1:k
        %    bar(t, c(:,j), 1, 'histc')
        %end
    end
  end
end

