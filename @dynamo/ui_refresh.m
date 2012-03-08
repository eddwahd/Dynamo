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

    if 0
        ax = subplot(2, 1, 1);
        self.plot_seq(ax, full);

        ax = subplot(2, 1, 2);
        self.plot_X(ax, full);
        %self.plot_stats(ax);
    else
        ax = get(h, 'CurrentAxes');
        self.plot_seq(ax, full);
    end
end
