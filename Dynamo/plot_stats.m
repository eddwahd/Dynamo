function plot_stats(stats, handle)
% Plots the statistics of an optimization run.

% Ville Bergholm 2011

pos = get(handle, 'position');

h = (pos(4)-pos(2))/2
dy = pos(4)/2; 

f = @semilogy; % @plot;

%iter = 1:length(stats.error);
%ax = plotyy(iter, abs(stats.error), iter, stats.wall_time, f);
% plotyy seems to have a bug related to the tickmarks of log plots

subplot('Position', [pos(1) 2*pos(2)+h pos(3) h]) 
f(stats.wall_time, abs(stats.error));
xlabel('Wall time (s)')
ylabel('Normalized error')
grid on

subplot('Position', [pos(1) pos(2) pos(3) h]) 
f(stats.wall_time, stats.fluence)
xlabel('Wall time (s)')
ylabel('Fluence')
grid on

%hold(ax(2), 'on');
%f(ax(2), iter, stats.cpu_time, 'r-');

%xlabel('Iteration')
%ylabel(ax(1), 'Error')
%ylabel(ax(2), 'Time (s)')
%legend(ax(2), 'Wall time', 'CPU time')
%a = axis(); a(4) = 1;
%axis(a)
