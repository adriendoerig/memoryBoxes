% temporal parameters
dt = 0.001;                 % [s]
boxLength = 0.2;            % [s]
tauIntegrate = 0.05;        % [s]
tauDecay = 1;               % [s]
readoutTime = 0.8;          % [s]
simulationTime = 1;         % [s]

% wongWang params
wongWang_gain = 1;          % gain from boxes stage to decision stage
wongWang_sigma = .01;       % noise from boxes stage to decision stage
wongWang_tStab = .5;       % [s], time to let the net stabilize
wongWang_mu0 = 25;          % wongWang "reactivity" -> high mu= = "jumpy" network

% a few commonstimuli to put in the boxes
vernier = ones(1,boxLength/dt);         
grating = zeros(1,boxLength/dt);
antivernier = -1*ones(1,boxLength/dt);

% declaring variables
stimulus = cell(1,nBoxes);  % one cell to describe each boxe's content
stimulus{1} = vernier;      % box 1 content
stimulus{2} = grating;      % box 2 content
stimulus{3} = antivernier;  % box 3 content

% compute dynamics of the memory boxes 
summedBoxOutputs = memoryBoxesDynamics( stimulus, tauIntegrate, tauDecay, readoutTime, simulationTime, dt);

% take boxes stage output and add gain and noise
wongWang_input = normrnd(wongWang_gain*summedBoxOutputs, wongWang_sigma);

% feed this to wongWang
[decision, DT, ~] = WongWangNew(wongWang_input, wongWang_tStab, wongWang_mu0);