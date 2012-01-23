function analyze()
% Analyzes the results of an optimization run.

% Ville Bergholm 2011


global OC;

err = OC.config.error_func();
fprintf('Final normalized error: %g\n    Wall time: %g s\n    CPU  time: %g s\nTermination reason: %s\n\n\n', ...
	err, OC.stats.wall_time(end), OC.stats.cpu_time(end), OC.opt.term_reason);

fprintf('Number of gradient evaluations: %d\n', OC.opt.N_eval);
fprintf('Final sequence duration: %g\n', sum(OC.seq.tau));


% rough error estimates
d = real(eig(OC.system.A));
n = length(d);
T = sum(OC.seq.tau);
e_max = 1-exp(-T*sum(d)/n)
e_min = 1-sum(exp(-T*d))/n


figure()
subplot(2, 1, 1)
plot_seq(OC.seq)

subplot(2, 1, 2)
plot_stats(OC.stats)
end
