function plot_seq(seq)
% Plots a control sequence using superposed bars.

% Ville Bergholm 2011


% start times for pulses
t = [0; cumsum(seq.tau)];

c = seq.control;
nc = size(c, 2); % number of controls

colormap(jet)
set(gca, 'CLim', [1 nc])
hold on
for j=1:length(t)-1
    x = [t(j), t(j+1), t(j+1), t(j)];
    
    [dummy, I] = sort(abs(c(j, :)), 2, 'descend'); % dummy instead of ~ for Octave
    for k=1:nc
        y = c(j, I(k));
        p = patch(x, [0, 0, y, y], 1);
        set(p, 'FaceColor','flat',  'CData',I(k), 'CDataMapping','scaled')
    end
end

axis tight
title('Control sequence')
xlabel('Time')
ylabel('Control amplitude')
grid on


%stairs(t, c)

%c = [c; zeros(1,nc)]; % final time step is a dummy%for j=1:k
%    bar(t, c(:,j), 1, 'histc')
%end
end
