function set_plotstyle(ax)
% Given an axes handle, sets the plotting style -related properties for those axes.
%
% If everything works properly, all the child elements of the
% axes (titles, labels, legends, tick labels etc.) should inherit
% these settings.
    
% Ville Bergholm 2012


set(ax, 'FontSize',18); %, 'FontName','Bitstream Vera Sans');
set(get(ax, 'Parent'), 'DefaultLineLineWidth',1); % apparently not an axes property(!)
set(ax, 'LineStyleOrder', '-|-.')
end
