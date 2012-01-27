function animate(self)
% Visualizes the evolution of the initial state under the current control sequence.

% TODO for now it only handles state ops in vec representation
    
    figure()
    s = state(inv_vec(self.X(0)));
    plot(s);
    a = axis(); a(end) = 1; % fix the z scaling 
    
    for k = 1:length(self.seq.tau)
        pause(0.1);
        k,
        s = state(inv_vec(self.X(k)));
        plot(s);
        axis(a); 
    end
end
