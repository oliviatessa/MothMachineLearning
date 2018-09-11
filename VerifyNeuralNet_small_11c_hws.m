%8/31/18
%Script to verify Callin's neural net predictions

%11c is for vertical aggressive maneuver

%% We the People, in Order to form a perfect Union,...
clc;
clearvars;
close all;
start_time = now;
disp(['Time at start: ',datestr(start_time)])
tic

%% Variable setup

LengthScaleFactor = 1; %Make sure this is an integer!
%Pack the parameter struct (NOT to be altered)
par.L1 = LengthScaleFactor*0.908; %Length from the thorax-petiole joint to 
    %the center of the head-thorax in cm - THIS WILL CHANGE 
          %FOR DAUBER (0.47195), BEE (0.2915), MOTH SCENARIO (0.908)
par.L2 = LengthScaleFactor*1.7475; %Length from the thorax-petiole joint 
    %to the center of the gaster in cm - THIS WILL CHANGE 
          %FOR DAUBER (1.177466667), BEE (0.2935), MOTH SCENARIO (1.7475)
par.L3 = LengthScaleFactor*0.75; %Length from the thorax-petiole joint to 
    %the aerodynamic force vector in cm - THIS WILL CHANGE 
          %FOR DAUBER (0.25), BEE (0.25), MOTH SCENARIO (0.75)
par.rho = 1; %The density of the insect in g/(cm^3)
par.rhoA = 1.18*10^-3; %The density of air in g/(cm^3)
par.muA = 1.86*10^-4; %The dynamic viscosity of air at 27C in g/(cm*s)
par.ahead = LengthScaleFactor*0.908; %Major axis of head-thorax ellipsoid 
          %in cm - THIS WILL CHANGE 
          %FOR DAUBER (0.47195), BEE (0.2915), MOTH SCENARIO (0.908)
par.abutt = LengthScaleFactor*1.7475; %Major axis of abdomen ellipsoid
          %in cm  - THIS WILL CHANGE 
          %FOR DAUBER (0.35), BEE (0.2935), MOTH SCENARIO (1.7475)
par.bhead = LengthScaleFactor*0.507; %Minor axis of head-thorax ellipsoid 
          %in cm - THIS WILL CHANGE 
          %FOR DAUBER (0.1839665), BEE (0.1815), MOTH SCENARIO (0.507)
par.bbutt = LengthScaleFactor*0.1295; %Minor axis of abdomen ellipsoid 
          %in cm - THIS WILL CHANGE 
          %FOR DAUBER (0.15), BEE (0.173), MOTH SCENARIO (0.1295)
par.K = LengthScaleFactor*29.3; %K is the torsional spring constant of 
          %the thorax-petiole joint in (cm^2)*g/(rad*(s^2))
par.c = LengthScaleFactor*14075.8; %c is the torsional damping constant of 
            %the thorax-petiole joint in (cm^2)*g/s
par.g = 980; %g is the acceleration due to gravity in cm/(s^2)

%% Import Callin's data set
NNpreds_small = csvread('CallinsData/NNpreds_small.csv',1);

%The 19 columns are in the following order:
%1 through 8
%x_0, y_0, phi_0, theta_0, x_dot_0, y_dot_0, phi_dot_0, theta_dot_0, 
%9 through 16
%F_pred, alpha_pred, tau_pred, x_99, y_99, phi_99, theta_99, x_dot_99_pred, 
%17 through 19
%y_dot_99_pred, phi_dot_99_pred, theta_dot_99_pred

ICs(:,1:2) = NNpreds_small(:,1:2); %x and y
ICs(:,3) = NNpreds_small(:,4); %theta
ICs(:,4) = NNpreds_small(:,3); %phi
ICs(:,5:6) = NNpreds_small(:,5:6); %x dot and y dot
ICs(:,7) = NNpreds_small(:,8); %theta dot
ICs(:,8) = NNpreds_small(:,7); %phi dot

%This is ordered x, y, theta, phi, xdot, ydot, theta dot, phi dot

F_inputs = NNpreds_small(:,9);
alpha_inputs = NNpreds_small(:,10);
tau_inputs = NNpreds_small(:,11);
%This is ordered F, alpha, tau

NumOfRuns = numel(NNpreds_small(:,1)); %This counts how many runs to do

%% More variables NOT to alter
hwbf = 50; %wing beat frequency of the insect in Hz - THIS WILL CHANGE 
          %FOR DAUBER (100), BEE (400), MOTH SCENARIO (50)
timestep = 100; %The number of timesteps we will use when interpolating
                %DO NOT CHANGE TIMESTEP
kk = 0; %Overall counter. Do not reset!
kkskip = 1; %Overall counter OF SKIPS. Do not reset!
nn = 1; %Counter of full runs. Do not reset!

%% Time vector
Tstore = linspace(0,((1/hwbf)*500),...
    (timestep*500))';
% save('11h-AMSumOfPrimes/Tstore_hws_sp.mat','Tstore');

timechunk = Tstore(100) - Tstore(1);

%% Prescribe the storage data
xx = zeros(timestep,8);
% Q = [];
tt = linspace(0,1/hwbf,timestep)';

lastQ = zeros(NumOfRuns,8);

timing1 = toc;

%% Using ode45 
tic

for i = 1:NumOfRuns
    
    % Initial conditions (MOVE INTO LOOP)
    %q0 = [x, y, theta, phi, xdot, ydot, thetadot, phidot]
    q0 = ICs(i,:);
    F = F_inputs(i,1);
    alpha = alpha_inputs(i,1);
    tau = tau_inputs(i,1);
    %Options for ode45
    OPTIONS = odeset('RelTol',1e-2,'AbsTol',1e-4); 

    par.betaR = q0(4) - q0(3) - pi; %This is the resting configuration of our 
    %torsional spring(s) = Initial abdomen angle - initial head angle - pi

    %To clear the pesky warnings/errors
    clc;

    %[TOUT,YOUT,TE,YE,IE] = ode45(ODEFUN,TSPAN,Y0,OPTIONS)
    [T, Q] = ode45(@(t,Q) myODE_4(t, Q, par, F, alpha,...
        tau), [0, 1/hwbf], q0, OPTIONS);
        
        xx = interp1(T, Q, tt);
        bigQ{i,1} = xx;
        lastQ(i,:) = xx(end,:);
        
disp(['Time at start: ',datestr(start_time)])
disp('Simulation in progress...')
disp(['Full run ', num2str(i), ' of ', num2str(NumOfRuns)]);
disp(['You are ',num2str(100*i/(NumOfRuns)), '% of the way there!']);

end 

%Save bigQ for personal verification
save('CallinsData/bigQ_small_hws_sp_con.mat','bigQ');

timing2 = toc;

%Create output file for Callin
NNpreds_output_small(:,1:2) = lastQ(:,1:2); %x and y
NNpreds_output_small(:,3) = lastQ(:,4); %phi
NNpreds_output_small(:,4) = lastQ(:,3); %theta
NNpreds_output_small(:,5:6) = lastQ(:,5:6); %x dot and y dot
NNpreds_output_small(:,7) = lastQ(:,8); %phi dot
NNpreds_output_small(:,8) = lastQ(:,7); %theta dot

csvwrite('CallinsData/NNpreds_small_outputJustNumbers.csv',NNpreds_output_small)

%% Write csv file
tic
csvHeader = {'x_0','y_0','phi_0','theta_0','x_dot_0','y_dot_0',...
    'phi_dot_0','theta_dot_0','F_pred',	'alpha_pred', 'tau_pred','x_99',...
    'y_99','phi_99','theta_99','x_dot_99_pred','y_dot_99_pred',...
    'phi_dot_99_pred','theta_dot_99_pred','x_a','y_a','phi_a','theta_a',...
    'x_dot_a','y_dot_a','phi_dot_a','theta_dot_a'};

commaHeader = [csvHeader;repmat({','},1,numel(csvHeader))]; %insert commaas
commaHeader = commaHeader(:)';
textHeader = cell2mat(commaHeader); %cHeader in text with commas

%write header to file
fid = fopen('CallinsData/NNpreds_small_output_all.csv','w'); 
fprintf(fid,'%s\n',textHeader);
fclose(fid);

NNpreds_output = [NNpreds_small,NNpreds_output_small];

%write data to end of file
dlmwrite('CallinsData/NNpreds_small_output_all.csv',NNpreds_output,...
    '-append');

timing3 = toc;

disp(['Loading vars took ', num2str(timing1), ' seconds'])
disp(['Running the verification took ',num2str(timing2),' seconds'])
disp(['Writing the csv file took ',num2str(timing3),' seconds'])

%% Final display on window
disp(' ');
disp('Completely done, son.');
disp(['Time at end: ',datestr(now)])