function ui_open(self)
% Create the UI figure.
    
    self.opt.UI_fig = figure('Name', 'Dynamo: Optimizing...', 'CloseRequestFcn', {@close_req, self});
    ax = axes();
    self.ui_refresh(true); % initial draw

    % FIXME
    drawnow();
    pause(3);
end


function close_req(src, event, self)
% Callback for closing the UI figure window.
    delete(src); % close the figure
    self.opt.UI_fig = [];
    self.opt.stop = true; % signal monitor_func that we should stop
end
