function analyze(self)
% Analyzes the results of an optimization run.


err = self.config.error_func(self);
fprintf('Final normalized error: %g\n    Wall time: %g s\n    CPU  time: %g s\nTermination reason: %s\n\n\n', ...
	err, self.stats.wall_time(end), self.stats.cpu_time(end), self.opt.term_reason);

fprintf('Number of gradient evaluations: %d\n', self.opt.N_eval);
fprintf('Final sequence duration: %g\n', sum(self.seq.tau));


% rough error estimates
d = real(eig(self.system.A));
n = length(d);
T = sum(self.seq.tau);
e_max = 1-exp(-T*sum(d)/n)
e_min = 1-sum(exp(-T*d))/n


figure()
subplot(2, 1, 1)
self.seq.plot()

subplot(2, 1, 2)

%handle = gca();
%pos = get(handle, 'position');
%h = (pos(4)-pos(2))/2;
%dy = pos(4)/2; 
%subplot('Position', [pos(1) 2*pos(2)+h pos(3) h]) 
%subplot('Position', [pos(1) pos(2) pos(3) h]) 

f = @loglog; %@semilogy % @plot;
[ax, h1, h2] = plotyy(self.stats.wall_time, abs(self.stats.error), self.stats.wall_time, self.stats.integral, f);
xlabel('Wall time (s)')
set(get(ax(1),'Ylabel'),'String','Normalized error') 
set(get(ax(2),'Ylabel'),'String','Control integral') 
grid on
set(h2,'LineStyle','--')


end
