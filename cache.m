classdef cache < matlab.mixin.Copyable
% Copyable handle class for doing the heavy computing and storing the results.

% Shai Machnes   2010-2011
% Ville Bergholm 2011-2012
    

  properties (SetAccess = private)
      H  % Generator for a time slice.
      P  % Propagator for a time slice. Roughly P = expm(-dt * H). Computed by calcPfromHfunc.
      U  % Forward propagators. U{k+1} = P{k} * U{k}
         % U{k} is the system at t = sum(tau(1:(k-1))) = t_{k-1}
      L  % Backward propagators. L{k-1} = L{k} * P{k-1};
         % L{k} is the adjoint system at t = sum(tau(1:(k-1))) = t_{k-1}
      g
      H_v  % eigendecomposition data for -dt*H, updated when P is updated
      H_eig_factor
  end

  properties (Access = public)
      H_needed_now  % flags
      P_needed_now
      U_needed_now
      L_needed_now
      g_needed_now
      
      E
      VUh
  end

  properties (Access = private)
      H_is_stale  % flags
      P_is_stale
      U_is_stale
      L_is_stale
      g_is_stale
      
      calcPfromHfunc
      UL_mixed  % flag: use mixed states in Hilbert space representation?
  end

  methods
      function self = cache(n_timeslots, U_start, L_end, use_eig, UL_hack)
      % Set up caching (once we know the number of time slices and U and L endpoints).

          temp = [1, n_timeslots];
          
          if use_eig
              % Store the eigendecomposition data as well
              self.H_v          = cell(temp);
              self.H_eig_factor = cell(temp);
              self.calcPfromHfunc = @calcPfromH_exact_gradient;
          else
              self.calcPfromHfunc = @calcPfromH_expm;
          end
          self.UL_mixed = UL_hack;

          self.H = cell(temp);
          self.P = cell(temp);
          self.U = cell(temp + [0, 1]);
          self.L = cell(temp + [0, 1]);

          % U: X_initial propagated forward up to a time instant.
          self.U{1} = U_start;
          % L: X_final' propagated backward, or in open-system tasks, a pure propagator (even though this requires more memory).
          self.L{end} = L_end;

          self.g = NaN;
          self.E = NaN;
          
          % Keep track of what needs re-computation if we want a complete update of everything.
          % U{1} and L{end} are never stale and never recomputed.
          self.H_is_stale = true(temp);
          self.P_is_stale = true(temp);
          self.U_is_stale = [false, true(temp)]; % Updates for H via 'control_update' get propagated automatically
          self.L_is_stale = [true(temp), false];
          self.g_is_stale = true;
          
          % Here we indicate which values we need to have up-to-date
          % The 'needed_now' function looks at everything we need to recompute now, 
          % compared to what in principle is_stale, and (in principle) executes the
          % optimal set of operations so that for everything which was marked
          % 'needed_now' is up-to-date

          % Example: We modified H{3}, so in theory we need to recompute U{4:end}.
          % But for our immediate needs we only want U{7} and L{7},
          % so we mark U{4:end} as "is_stale", but only 4:7 as "needed_now".

          self.H_needed_now = false(temp);   
          self.P_needed_now = false(temp);   
          self.U_needed_now = [false, false(temp)];
          self.L_needed_now = [false(temp), false]; 
          self.g_needed_now = false;
      end
      
      
      function invalidate(self)
      % Invalidates the entire cache.

          self.H_is_stale(:) = true;
          self.P_is_stale(:) = true;
          self.U_is_stale(2:end) = true;
          self.L_is_stale(1:(end-1)) = true;
          self.g_is_stale = true;
      end


      function mark_as_stale(self, changed_t_mask)
      % Marks the selected timeslots as stale.

          self.H_is_stale(changed_t_mask) = true;
          self.P_is_stale(changed_t_mask) = true;

          % Propagate the H_is_stale to the U and Ls.
          self.U_is_stale( (find(self.H_is_stale, 1, 'first')+1):end) = true;
          self.L_is_stale(1:find(self.H_is_stale, 1, 'last'))         = true;
          self.g_is_stale = true;
      end
      
      
      function refresh(self, sys, tau, fields)
      % This function does most of the heavy computing.
      % It should be called _after_ setting _all_ the required *_needed_now fields
      % to minimize unnecessary computational work.
      %
      % The function looks at everything we need to recompute now,
      % compared to what in principle is_stale, and (in principle) executes the
      % optimal set of operations so that for everything which was marked
      % 'needed_now' is up-to-date.
      %
      % It then updates the 'is_stale' to indicate what has been actually
      % computed, and clears the 'needed_now'
      %
      % The inter-dependence of H and U/L updates is taken care of in mark_as_stale()

          n_timeslots = length(self.H);
          
          % computing g may require additional U and L elements, so check that first
          g_recompute_now = self.g_needed_now && self.g_is_stale;
          if g_recompute_now
              % g can be computed using any slice k \in [1, n+1]: g = trace(L_k * U_k).
              % Try to figure out which k requires least additional computation.
              g_k = self.g_setup_recalc();
          end
          
          U_recompute_now = self.U_needed_now & self.U_is_stale;
          L_recompute_now = self.L_needed_now & self.L_is_stale;
          % To recompute U, you need to start at a cell that is fully recomputed.
          % U{1} and L{n_timeslots+1} are by definition never stale and never recomputed.
          for t=n_timeslots:(-1):2
              if U_recompute_now(t+1) && self.U_is_stale(t)
                  U_recompute_now(t) = true;
              end
          end
          for t=2:n_timeslots
              if L_recompute_now(t-1) && self.L_is_stale(t)
                  L_recompute_now(t) = true;
              end
          end

          % Now that we know which Us & Ls we need to recompute, we can figure out which Ps and Hs must be up-to-date
          P_recompute_now = (U_recompute_now(2:end) | L_recompute_now(1:(end-1)) | self.P_needed_now) & self.P_is_stale;
          H_recompute_now = (P_recompute_now | self.H_needed_now) & self.H_is_stale;

          % Compute the Hamiltonians
          h_idx = find(H_recompute_now);
          for t=h_idx
              H = sys.A;
              for c = 1:length(sys.B)
                  u = fields(t, c);
                  H = H + u * sys.B{c};
              end
              self.H{t} = H;
          end

          % Compute the exp(H) and any other per-H computation which may be needed for the gradient function
          p_idx = find(P_recompute_now);
          for t=p_idx
              self.calcPfromHfunc(self, t, tau(t)); % Compute the Ps - a single piece of propagator
              % NOTE: calcPfromHfunc may also compute other values which will be needed for gradient calculations.
              %       These should be stored in cache. Their up-to-date-ness is identical to that of P.
          end

          % Compute the Us - forward propagation (we never recompute U{1})
          u_idx = find(U_recompute_now);
          % Compute the Ls - adjoint system propagation
          el_idx = fliplr(find(L_recompute_now));
          if self.UL_mixed
              % mixed states, unitary evolution: propagate from both sides
              for t=u_idx
                  self.U{t} = self.P{t-1} * self.U{t-1} * self.P{t-1}';
              end
              for t=el_idx
                  self.L{t} = self.P{t}' * self.L{t+1} * self.P{t};
              end
          else
              % propagate U from left, L from right
              for t=u_idx
                  self.U{t} = self.P{t-1} * self.U{t-1};
              end
              for t=el_idx
                  self.L{t} = self.L{t+1} * self.P{t};
              end
          end

          % and finally g
          if g_recompute_now
              self.g = trace_matmul(self.L{g_k}, self.U{g_k});
              self.g_is_stale = false;
          end
              
          % Mark what has been actually computed
          temp = [1, n_timeslots];

          self.H_is_stale(H_recompute_now) = false;
          self.P_is_stale(P_recompute_now) = false;
          self.U_is_stale(U_recompute_now) = false;
          self.L_is_stale(L_recompute_now) = false;
          
          self.H_needed_now = false(temp);
          self.P_needed_now = false(temp);
          self.U_needed_now = [false, false(temp)];
          self.L_needed_now = [false(temp), false];
          self.g_needed_now = false;

          %ret = [sum(H_recompute_now), sum(P_recompute_now), sum(U_recompute_now), sum(L_recompute_now)]
      end
  
      
      function calcPfromH_exact_gradient(self, t, dt)
      % Computes cache.P{t} using the eigendecomposition, stores some
      % extra stuff for cheap exact gradient computation later on.

          minus_dt_H = -dt * self.H{t};
          N = length(minus_dt_H);

          %% Compute the eigenvalue factors and eigenvectors of -dt*H

          [v, zeta, exp_d] = eig_factors(minus_dt_H, true);
          self.H_v{t} = v;
          self.H_eig_factor{t} = zeta;

          %% And finally expm(-dt*H) using the eigendecomposition

          self.P{t} = v * diag(exp_d) * v';
      end


      function calcPfromH_expm(self, t, dt)
      % Compute P{t} using expm.
          self.P{t} = expm(-dt * self.H{t});
      end
      
      
      function t = g_setup_recalc(self)
      % Returns the optimal time slice in which to compute g.

      % Future work: Optimally, we can search for the time which will require minimum calculations to compute the value.
      % However, this is tricky (we need to count the expensive H-->P computations (expm or similar) and maybe weigh-in the cheap
      % U and L updates (matrix multiplications).

      % But for now, we'll do a sub-optimal algorithm and take the left-most L that is up-to-date (and therefore update all L-s to the
      % right of it). If all the slots updated are in a single block, we'll be optimal

          flags = (~self.U_is_stale) + (~self.L_is_stale)*2;

          flag3 = find(flags==3,1,'first');
          if ~isempty(flag3)
              t = flag3;
          else
              p1 = find(flags>0);
              p2 = find((flags(p1(1:(end-1)))==1) & (flags(p1(2:end))==2));

              if isempty(p2)
                  error('There should be at least one');
              end
              cost = p1(p2+1)-p1(p2);
              %    [~,mincostpos] = min(cost);
              [dummy, mincostpos] = min(cost); % for Octave
              t = p1(p2(mincostpos));
          end
          self.U_needed_now(t) = true;
          self.L_needed_now(t) = true;
      end
  end
end
