function plot_stats(stats)
% Plots the statistics of an optimization run.

% Ville Bergholm 2011

iter = 1:length(stats.Q_func);

f = @semilogy; % @plot;

%ax = plotyy(iter, abs(1-stats.Q_func), iter, stats.wall_time, f);
% plotyy seems to have a bug related to the tickmarks of log plots

f(stats.wall_time, abs(1-stats.Q_func));
xlabel('Wall time (s)')
ylabel('Error')

%hold(ax(2), 'on');
%f(ax(2), iter, stats.cpu_time, 'r-');

%xlabel('Iteration')
%ylabel(ax(1), 'Error')
%ylabel(ax(2), 'Time (s)')
%legend(ax(2), 'Wall time', 'CPU time')
grid on
a = axis(); a(4) = 1;
axis(a)
title('Optimization statistics')
