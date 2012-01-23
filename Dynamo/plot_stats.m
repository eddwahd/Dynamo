function plot_stats(stats)
% Plots the statistics of an optimization run.

% Ville Bergholm 2011


%handle = gca();
%pos = get(handle, 'position');
%h = (pos(4)-pos(2))/2;
%dy = pos(4)/2; 
%subplot('Position', [pos(1) 2*pos(2)+h pos(3) h]) 
%subplot('Position', [pos(1) pos(2) pos(3) h]) 

f = @loglog; %@semilogy % @plot;
[ax, h1, h2] = plotyy(stats.wall_time, abs(stats.error), stats.wall_time, stats.fluence, f);
xlabel('Wall time (s)')
set(get(ax(1),'Ylabel'),'String','Normalized error') 
set(get(ax(2),'Ylabel'),'String','Control fluence') 
grid on
set(h2,'LineStyle','--')
