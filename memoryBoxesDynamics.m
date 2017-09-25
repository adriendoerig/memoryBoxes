function [ summedBoxOutputs ] = memoryBoxesDynamics( stimulus, tauIntegrate, tauDecay, readoutTime, simulationTime, dt)
%MEMORYBOXESDYNAMICS Manages the memory boxes and returns their output
% Parameters:
% inputSequence : vector of 1, 0 & -1: 1 for vernier, -1 for anti-vernier
% inputGain     : "strength" of input
% tau           : time constant
% dt            : time step [s]
% varargin      : you may include an initial  value (default = 0)% to play around

plottingMemoryTraces = 1;                           % plots box traces if set to 1.

nBoxes = length(stimulus);
boxLength = length(stimulus{1})*dt;                 % [s]
memoryTraces = zeros(nBoxes, simulationTime/dt);    % will store memory traces in each box
summedBoxOutputs = 0;                               % will store final result

% all boxes go
for i = 1:nBoxes
    % integration
    [integrationOutput,memoryTraces(i,boxLength/dt*(i-1)+1:boxLength/dt*i)] = boxIntegrate(stimulus{i},tauIntegrate,dt);

    % decay
    memoryTraces(i,boxLength/dt*i+1:end) = boxDecay(integrationOutput,simulationTime-i*boxLength, tauDecay, dt);
    
    % summing box states at readout time.
    summedBoxOutputs = summedBoxOutputs+memoryTraces(i,readoutTime/dt);
    
end

if plottingMemoryTraces
    figure(1000)
    for i = 1:nBoxes
        subplot(nBoxes,1,i)
        plot(1:simulationTime/dt, memoryTraces(i,:))
        line([readoutTime/dt, readoutTime/dt], [-.05 .05], 'Color', 'g');
        mtit(['Memory boxes in temporal order. Summed value for decision stage = ', num2str(summedBoxOutputs)])
    end
end

