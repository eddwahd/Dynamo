function analyze()
% Analyzes the results of an optimization run.

% Ville Bergholm 2011


global OC;


Q = OC.config.Q_func();
 
fprintf('Fidelity reached: 1 - %g\n    Wall time: %g s\n    CPU  time: %g s\nTermination reason: %s\n\n\n', ...
	1-Q, OC.stats.wall_time(end), OC.stats.cpu_time(end), OC.opt.term_reason);


% Compare normalized errors
%X = system_get(length(OC.seq.tau));  % X(t_n)

%d2 = norm(OC.system.X_final -X, 'fro')^2 / OC.system.norm2;
%B2 = norm(X, 'fro')^2 / OC.system.norm2;

%fprintf('|X_f|^2 = %g\nNormalized error squared = %e\nNormalized |X(t)|^2 = %g\n', ...
%         OC.system.norm2, d2, B2);

fprintf('Number of gradient evaluations: %d\n', OC.opt.N_eval);
fprintf('Final sequence duration: %g\n', sum(OC.seq.tau));


% rough error estimates
d = real(eig(OC.system.A));
n = length(d);
T = sum(OC.seq.tau);
e_max = 1-exp(-T*sum(d)/n)
e_min = 1-sum(exp(-T*d))/n


figure()
subplot(3, 1, 1)
plot_seq(OC.seq)
subplot(3, 1, 2)
plot_stats(OC.stats)
subplot(3, 1, 3)
plot(OC.stats.wall_time, OC.stats.fluence)

end

