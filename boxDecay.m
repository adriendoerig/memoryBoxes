function [ boxDecayOutput ] = boxDecay( boxOutput, decayTime, tau, dt )
%BOXINTEGRATE Decays memory box content after it has closed
% Parameters:
% boxOutput     : box content before decay
% decayTime     : duration of decay [s]
% tau           : time constant
% dt            : time step [s]

plotting = 0;   % set to 1 to plot output

boxDecayOutput = zeros(1,int32(decayTime/dt));

boxDecayOutput(1) = boxOutput;

for t = 1:decayTime/dt-1
    boxDecayOutput(t+1) = boxDecayOutput(t)*(1-dt/tau);
end

if plotting
    figure()
    plot(1:length(boxDecayOutput),boxDecayOutput)
    xlabel('time [s*dt]')
    ylabel('memory box content: decay')
end
    

end

