function [decision, DT, success] = WongWangNew(v, t_stab, mu0) %%% Wong & Wang, JNS 2006 appendix %%%
% decision = 1 for first vernier
% DT = decision time
% success = 1 if the network successfully reaches a decision.

% in case the network fucks up and doesn't reach a decision, 
% these are just average values so it doesn't bug during fminsearch.
decision = 0;
DT = 300;
success = 0;

plotting = 1; % plots network state when set to 1

%%%% Synaptic time and other constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tnmda = 100;   % NMDA decay time constant; JNS 2006
Tampa = 2;      % AMPA decay time constant
gamma = 0.641;  % Parameter that relates presynaptic input firing rate to synaptic gating variable
JAext = 0.00052; % Synaptic coupling constant to external inputs; JNS 2006

%%%% FI curve parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = 270; b = 108; d = 0.1540;  % Parameters for excitatory cells; see JNS2006

%%%% Vectorizing variables for evaluations at end of (block) loop %%%%%%%%%%

r1_traj = [];  r2_traj = []; % firing rates
s1_traj = [];  s2_traj = []; % synaptic gating variables

%%%% Important parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

threshold = 15;       % threshold [Hz]
% mu0       = 20.0;      % Stimulus strength = input received when there is no bias (LDD)
noise_amp = 0.02;      % Noise amplitude into selective populations
N_trials  = 1 ;       % Total number of trials

%%% Time conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
dt         = 0.5;                  % Time step (ms)
time_wind  = 2/dt;                 % Temporal window size for averaging
T_total    = 5000/dt + time_wind;  % Total number of steps
slide_wind = time_wind;             % Sliding step for window



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Trials number and (block) loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ww = 1:N_trials % Trial loop

    trial_no = ww; 

    %---- Vectorise variables with sliding window -------------
    
    nu1_wind = [] ; nu2_wind = [] ;
    s1_wind  = [] ; s2_wind  = [] ;

    %---- Initial conditions and clearing variables -----------
    
    s1_in=0.1; s2_in=0.1; 
    clear nu1_in nu2_in I_eta1 I_eta2 ;
    nu1_in = 2; nu2_in = 2; 
    I_eta1_in = noise_amp*randn ; I_eta2_in = noise_amp*randn ;

    %---- Intialise and vectorise variables to be used in loops below ------
    
    s1 = s1_in.*ones(1,T_total); s2 = s2_in.*ones(1,T_total);
    nu1 = nu1_in.*ones(1,T_total); nu2 = nu2_in.*ones(1,T_total); %firing rates?? (LDD)
    phi1 = nu1_in.*ones(1,T_total); phi2 = nu2_in.*ones(1,T_total);
    I_eta1 = I_eta1_in.*ones(1,T_total); I_eta2 = I_eta2_in.*ones(1,T_total);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Time loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

           
    for t = 1:T_total

        %---- Constant effective external current input (with inhibition taken into account)
        
        I0E1 = 0.3255; I0E2 = 0.3255;

        %---- Random dot stimulus------------------------------------------------------
        
        if t > t_stab/dt*1000
            I_stim_1 = (JAext*mu0*(1+v)); % To population 1
            I_stim_2 = (JAext*mu0*(1-v)); % To population 2
        else
            I_stim_1 = (JAext*mu0); % To population 1
            I_stim_2 = (JAext*mu0); % To population 2 
        end
        
        %---- Recurrent synaptic coupling constants-------------------------------------
        
        JN11 = 0.2609; JN22 = 0.2609; 
        JN12 = 0.0497; JN21 = 0.0497; 
        
        %---- Resonse function of competiting excitatory population 1 ------
        
        % Total synaptic input to population 1
        Isyn1(t) = JN11.*s1(t) - JN12.*s2(t) + I_stim_1 + I_eta1(t);
        
        % Transfer function of population 1
        phi1(t)  = (a.*Isyn1(t)-b)./(1-exp(-d.*(a.*Isyn1(t)-b)));

	    %---- Response function of competiting excitatory population 2 -----
	    
        % Total synaptic input to population 2
        Isyn2(t) = JN22.*s2(t) - JN21.*s1(t) + I_stim_2 + I_eta2(t);

        % Transfer function of population 2
	    phi2(t)  = (a.*Isyn2(t)-b)./(1-exp(-d.*(a.*Isyn2(t)-b)));

	    %---- Dynamical equations -------------------------------------------

	    % Mean NMDA-mediated synaptic dynamics updating
	    s1(t+1) = s1(t) + dt*(-s1(t)/Tnmda + (1-s1(t))*gamma*nu1(t)/1000);
	    s2(t+1) = s2(t) + dt*(-s2(t)/Tnmda + (1-s2(t))*gamma*nu2(t)/1000);

        % Ornstein-Uhlenbeck generation of noise in pop1 and 2
        I_eta1(t+1) = I_eta1(t) + (dt/Tampa)*(I0E1-I_eta1(t)) + sqrt(dt/Tampa)*noise_amp*randn ;
        I_eta2(t+1) = I_eta2(t) + (dt/Tampa)*(I0E2-I_eta2(t)) + sqrt(dt/Tampa)*noise_amp*randn ;

        % To ensure firing rates are always positive. Large noise amplitude 
        % may result in unwanted negative values
        if phi1(t) < 0 
            nu1(t+1) = 0;
            phi1(t) = 0;
        else
            nu1(t+1) = phi1(t);
        end;
        if phi2(t) < 0
            nu2(t+1) = 0;
            phi2(t) = 0;
        else
            nu2(t+1) = phi2(t);
        end;
        
	    %==============================================================================================
        
    end;  %---- End of time loop --------

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%    %---- Calculating the mean rates and gating variables with sliding window -----

    nu1_wind = [nu1_wind (mean(nu1(1:time_wind)))];
    nu2_wind = [nu2_wind (mean(nu2(1:time_wind)))];
    s1_wind  = [s1_wind (mean(s1(1:time_wind)))];
    s2_wind  = [s2_wind (mean(s2(1:time_wind)))];

    
    for jj=1:((T_total-time_wind)/slide_wind)
        nu1_wind = [nu1_wind mean(nu1(jj*slide_wind:(jj*slide_wind+time_wind)))];
        nu2_wind = [nu2_wind mean(nu2(jj*slide_wind:(jj*slide_wind+time_wind)))];
        s1_wind = [s1_wind mean(s1(jj*slide_wind:(jj*slide_wind+time_wind)))];
        s2_wind = [s2_wind mean(s2(jj*slide_wind:(jj*slide_wind+time_wind)))];
    end;

    nu1_wind = [mean(nu1(1:time_wind)) nu1_wind];
    nu2_wind = [mean(nu2(1:time_wind)) nu2_wind];

    r1_traj = [r1_traj; nu1_wind]; r2_traj = [r2_traj; nu2_wind];
    s1_traj = [s1_traj; s1_wind];  s2_traj = [s2_traj; s2_wind];

    clear nu1 nu2 s1 s2 nu1_wind nu2_wind s1_wind s2_wind ;
    clear phi1 phi2 I_eta1 I_eta2 I_eta1_in I_eta2_in ;


end; %---- End trial loop ---------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Plots of first trajectories %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
if plotting
    for ww=1:1
        figure('Position',[900 439 612 439])
        plot([(dt*time_wind-dt*slide_wind):dt*slide_wind:(dt*T_total-dt*time_wind)],r1_traj(ww,1:end-1));
        hold on;
        plot([(dt*time_wind-dt*slide_wind):dt*slide_wind:(dt*T_total-dt*time_wind)],r2_traj(ww,1:end-1),'r');
        hold on;
        plot([(dt*time_wind-dt*slide_wind):dt*slide_wind:(dt*T_total-dt*time_wind)],threshold*ones(length(r2_traj)-1),'g');
    end;
    hold off
    xlabel('Time (ms)');ylabel('Firing rate (Hz)'); grid on;
    legend('1st vernier', '2nd vernier')
end

for t = 1:length(r1_traj)
    if r1_traj(1,t)>=threshold
        decision = 1;
        DT = t*2 - t_stab*1000; % we want ms but there is one number in r1_traj every TWO millieconds
        success = 1; % just to check that a decision has been made
        break
    end
    if r2_traj(1,t)>=threshold
        decision = -1;
        DT = t*2 - t_stab*1000;
        success = 1;
        break
    end
end
