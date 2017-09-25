% to play around

plottingMemoryTraces = 1;   % plots box traces if set to 1.

% parameters
dt = 0.001;                 % [s]
boxLength = 0.2;            % [s]
tauIntegrate = 0.05;        % [s]
tauDecay = 1;               % [s]
readoutTime = 0.8;          % [s]
simulationTime = 1;

% a few stimuli to put in the boxes
vernier = ones(1,boxLength/dt);         
grating = zeros(1,boxLength/dt);
antivernier = -1*ones(1,boxLength/dt);

% decalring variables
nBoxes = 3;                         % number of independant boxes
stimulus = cell(1,nBoxes);          % one cell to describe each boxe's content
boxOutputs = cell(1,nBoxes);        % will store integration output for each box
boxDecayOutputs = cell(1,nBoxes);   % will store decaying stored info in each box
summedBoxOutputs = 0;               % this will be fed to wongWangNew
stimulus{1} = vernier;
stimulus{2} = grating; 
stimulus{3} = antivernier; 

% will store memory traces in each box
memoryTraces = zeros(nBoxes, simulationTime/dt);


% all boxes go
for i = 1:nBoxes
    i
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

% wongWang_input = normrnd(wongWang_gain*, sigmaW);
% 
% [decision, DT, ~] = WongWangNew(wongWang_input, t_stab, wongWang_mu0)