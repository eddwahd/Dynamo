function ui_refresh(self, full, err)
% Refreshes the UI figure, plots stuff.
    
    % Do a full refresh or just redraw the graphical objects?
    if nargin < 2
        full = true;
    end
    
    h = self.opt.ui_fig;
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
        title(ax, self.system.description);
    end
    if nargin == 3
        % TODO Frobenius error
        %err = sqrt(2 * err * self.system.norm2);
        text(0.05, 0.9, sprintf('Error: %6.6g', err), 'Units','normalized',...
             'FontSize',18, 'BackgroundColor',[.8 .8 1])
    end

    drawnow();
end
