%% Demo 1 for the DYNAMO Quantum Optimal Control platform
%
% DYNAMO - Quantum Dynamic Optimization Package
% (c) Shai Machnes 2010, Institute of Theoretical Physics, Ulm University, Germany
% email: shai.machnes at uni-ulm.de
%
% All computer programs/code/scripts are released under the terms of the 
% GNU Lesser General Public License 3.0, except where explicitly 
% stated otherwise. Everything else (documentation, algorithms, etc) is 
% licensed under the Creative Commons Attribution-Share Alike 3.0 License
% (see "LICENSE.txt" for details). 
%
% For the latest version, visit http://www.qlib.info
%

%% What this demo does

clc

fprintf ('This demo will optimize a simple two-qubit QFT gate generation problem using a wide variety of algorithms\n');
fprintf ('\n');
fprintf ('If your interests are focued on finding optimal control sequences for your specific system, please refer to demo 00.\n');
fprintf ('However, if you are interested in OC algorithm research, this is the place for you.\n');
fprintf ('\n\n');


%% Now some search methods


% Which time slots do you want to modify ?
controls_mask = dyn.full_mask();
controls_mask(40:45, :) = false; % Keep these at initial random value, for no good reason other than demoing capabilities

% Try various methods ------------------------------------------------------------------

termination_conditions = struct( ...
    'loop_count',           1e10, ...
    'goal',                 1 - 1e-6, ...
    'wall_time_to_stop',    60, ...
    'cputime_to_stop',      60, ...
    'gradient_norm_min',    1e-20);

wall0 = now(); cpu0 = cputime(); fprintf ('Krotov (1st order update scheme, serial timeslot update - 1 step per timeslot)\n'); drawnow;
termination_reason = Krotov_search_function (controls_mask, termination_conditions);
fprintf('Goal reached: 1 - %g\n    Wall time: %g\n    CPU time:  %g\nTermination reason: %s\n\n\n', 1-get_current_value(),(now()-wall0)*(24*60*60), cputime()-cpu0, OC.const.termination_reason_str{termination_reason});
invalidate_timeslots_calc_cache(); intialize_timeslot_controls (initial_controls);

wall0 = now(); cpu0 = cputime(); fprintf ('1st order update scheme, modifying all timeslices concurrently (at each step)\n'); drawnow;
termination_reason = First_order_search_function (controls_mask,termination_conditions);
fprintf('Goal reached: 1 - %g\n    Wall time: %g\n    CPU time:  %g\nTermination reason: %s\n\n\n', 1-get_current_value(),(now()-wall0)*(24*60*60), cputime()-cpu0, OC.const.termination_reason_str{termination_reason});
invalidate_timeslots_calc_cache(); intialize_timeslot_controls (initial_controls);

wall0 = now(); cpu0 = cputime(); fprintf ('GRAPE (BFGS, all timeslices)\n'); drawnow;
OC.config.BFGS = struct('fminopt', struct('Display', 'off'));
termination_reason = BFGS_search_function (controls_mask,termination_conditions);
fprintf('Goal reached: 1 - %g\n    Wall time: %g\n    CPU time:  %g\nTermination reason: %s\n\n\n', 1-get_current_value(),(now()-wall0)*(24*60*60), cputime()-cpu0, OC.const.termination_reason_str{termination_reason});
invalidate_timeslots_calc_cache(); intialize_timeslot_controls (initial_controls);

% And now some hybrid schemes

per_block_termination_conditions = struct( ...
    'loop_count',           10, ...
    'goal',                 1 - 1e-6, ...
    'wall_time_to_stop',    6, ...
    'cputime_to_stop',      6, ...
    'gradient_norm_min',    1e-20);
meta_block_termination_conditions = struct( ...
    'loop_count',           1e4, ...
    'goal',                 1 - 1e-6, ...
    'wall_time_to_stop',    60, ...
    'cputime_to_stop',      60, ...
    'gradient_norm_min',    1e-20);

wall0 = now(); cpu0 = cputime(); fprintf ('Hybrid scheme: 1st order update method, updating 12 timeslices on each step, 10 steps per block before moving on to the next 12 timeslices\n'); drawnow;
termination_reason = Block_cycle_search_function (controls_mask, meta_block_termination_conditions, per_block_termination_conditions, @First_order_search_function, 12);
fprintf('Goal reached: 1 - %g\n    Wall time: %g\n    CPU time:  %g\nTermination reason: %s\n\n\n', 1-get_current_value(),(now()-wall0)*(24*60*60), cputime()-cpu0, OC.const.termination_reason_str{termination_reason});
invalidate_timeslots_calc_cache(); intialize_timeslot_controls (initial_controls);

first_stage_termination_condition = termination_conditions;
first_stage_termination_condition.goal = 1 - 0.1;
wall0 = now(); cpu0 = cputime(); fprintf ('Crossover demo: start with Krotov, finish with GRAPE. Crossover on goal of 1-%g\n', 1-first_stage_termination_condition.goal); drawnow;
termination_reason = Two_method_crossover_function (controls_mask, first_stage_termination_condition, @Krotov_search_function, termination_conditions, @BFGS_search_function);
fprintf('Goal reached: 1 - %g\n    Wall time: %g\n    CPU time:  %g\nTermination reason: %s\n\n\n', 1-get_current_value(),(now()-wall0)*(24*60*60), cputime()-cpu0, OC.const.termination_reason_str{termination_reason});
invalidate_timeslots_calc_cache(); intialize_timeslot_controls (initial_controls);

fprintf('\n\nAll done\n');
tmp=load('handel.mat'); sound(tmp.y(1600:18000), tmp.Fs);
