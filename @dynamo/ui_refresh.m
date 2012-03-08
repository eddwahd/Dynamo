function ui_refresh(self, full)
% Refreshes the UI figure, plots stuff.
    
    % Do a full refresh or just redraw the graphical objects?
    if nargin < 2
        full = true;
    end
    
    h = self.opt.UI_fig;
    % It's incredible how much work it takes just to make
    % MATLAB not steal window focus when it plots something.
    set(0, 'CurrentFigure', h);
    %ax = get(h, 'CurrentAxes');
    %self.seq.plot(ax, true);

    ax = subplot(2, 1, 1);
    self.plot_seq(ax, full);

    ax = subplot(2, 1, 2);
    self.plot_X(ax);
    %self.plot_stats(ax);
end
