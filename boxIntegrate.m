function [ boxOutputValue, boxOutputTrace ] = boxIntegrate( inputSequence, tau, dt, varargin )
%BOXINTEGRATE Integrates input in a memory box
% Parameters:
% inputSequence : vector of 1, 0 & -1: 1 for vernier, -1 for anti-vernier
% inputGain     : "strength" of input
% tau           : time constant
% dt            : time step [s]
% varargin      : you may include an initial  value (default = 0)

plotting = 0;   % set to 1 to plot output

boxOutputTrace = zeros(1,length(inputSequence));
    
if nargin == 4
    boxOutputTrace(1) = varargin{1};
else
    boxOutputTrace(1) = 0;
end

for t = 1:length(inputSequence)-1
    boxOutputTrace(t+1) = boxOutputTrace(t)*(1-dt/tau) + inputSequence(t)*dt;
end

boxOutputValue = boxOutputTrace(end);

if plotting
    figure()
    plot(1:length(inputSequence),boxOutputTrace)
    xlabel('time [s*dt]')
    ylabel('memory box content: integration')
end
    

end

