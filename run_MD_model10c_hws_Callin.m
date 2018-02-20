%4/4/17 Setting up my model (Model 1) in MATLAB using ode45
%Many models since then (and even more errors)
%Model 8 was essentially refined on 6/29/17
%7/6/17- Model 9 development will incorporate our cost function and 
%multiple wing strokes.
%8/7/17- Model 10(a?) will do all of the above for SEVERAL full run 
%throughs. Model 10(b?) will have a Gaussian cost function. a.k.a. no more
%window!
%9/4/17- FIXED K and c
%9/6/17- Fixed tau to have a maximum value of 10000 g*(cm^2)/(s^2)

%10b is for horizontal flight
%10c is for aggressive maneuver
%10d is for hovering
%10e is for horizontal flight, with restricted F and alphas after 5th ws
%10f is for aggressive maneuver, with restricted F and alphas after 5th ws
%10g is for hovering, with restricted F and alphas after 5th ws
%10h is for sum of prime number sines

%% We the People, in Order to form a perfect Union,...
clc;
clear all;
close all;

%% Variables

global L1 L2 L3 rho rhoA muA ahead abutt bhead bbutt K c F alpha tau0 f...
    g thetaR

L1 = 0.908; %Length from the thorax-petiole joint to the center of the 
    %head-thorax in cm - THIS WILL CHANGE 
          %FOR DAUBER (0.47195), BEE (0.2915), MOTH SCENARIO (0.908)
L2 = 1.7475; %Length from the thorax-petiole joint to the center of the 
    %gaster in cm - THIS WILL CHANGE 
          %FOR DAUBER (1.177466667), BEE (0.2935), MOTH SCENARIO (1.7475)
L3 = 0.75; %Length from the thorax-petiole joint to the aerodynamic force 
    %vector in cm - THIS WILL CHANGE 
          %FOR DAUBER (0.25), BEE (0.25), MOTH SCENARIO (0.75)
rho = 1; %The density of the insect in g/(cm^3)
rhoA = 1.18*10^-3; %The density of air in g/(cm^3)
muA = 1.86*10^-4; %The dynamic viscosity of air at 27C in g/(cm*s)
ahead = 0.908; %For Dauber %Major axis of head-thorax ellipsoid 
          %in cm - THIS WILL CHANGE 
          %FOR DAUBER (0.47195), BEE (0.2915), MOTH SCENARIO (0.908)
abutt = 1.7475; %Major axis of gaster ellipsoid in cm  - THIS WILL CHANGE 
          %FOR DAUBER (0.35), BEE (0.2935), MOTH SCENARIO (1.7475)
bhead = 0.507; %Minor axis of head - thorax ellipsoid in cm - THIS WILL CHANGE 
          %FOR DAUBER (0.1839665), BEE (0.1815), MOTH SCENARIO (0.507)
bbutt = 0.1295; %Minor axis of gaster ellipsoid in cm - THIS WILL CHANGE 
          %FOR DAUBER (0.15), BEE (0.173), MOTH SCENARIO (0.1295)
K = 29.3; %K is the torsional spring constant of the thorax-petiole joint
          %in (cm^2)*g/(rad*(s^2))
c = 14075.8; %c is the torsional damping constant of the thorax-petiole joint
             %in (cm^2)*g/s
% F = 243; %243; %F is the aerodynamic vector in cm/(s^2)
% alpha = pi/4; %alpha is the angle of the aerodynamic vectory with 
% %respect to the head-thorax
% tau0 = 0; %The initial torque applied
f = 4; %frequency in Hz
%tau = tau0*sin(2*pi*f*tspan); %tau is the torque applied
g = 980; %g is the acceleration due to gravity in cm/(s^2)

shift = 0; %This is to increase the file number as appropriate

hwbf = 50; %wing beat frequency of the insect in Hz - THIS WILL CHANGE 
          %FOR DAUBER (50), BEE (200), MOTH SCENARIO (25)
% flyspd = 300; %Flying speed of the insect in cm/s - THIS WILL CHANGE
%           %FOR DAUBER (?), BEE (670.56), MOTH (300).
% bl = 5.311; %The body length of the organism in cm - THIS WILL CHANGE
%           %FOR DAUBER (2.31), BEE (1.153), MOTH (5.311).
timestep = 100; %The number of timesteps we will use when interpolating
                %DO NOT CHANGE TIMESTEP
halfwingStrokes = 1; 
% numOfTrajectories = 1000; %Number of trajectories for a given wing beat
FullRuns = 3; %Number of full runs 
kk = 1; %Overall counter. Do not reset!
kkskip = 1; %Overall counter OF SKIPS. Do not reset!
nn = 1; %Counter of full runs. Do not reset!

Tstore = zeros(timestep*halfwingStrokes, 1);
Tstore = linspace(0,((1/hwbf)*halfwingStrokes),...
    (timestep*halfwingStrokes))';
% save('1c-AggressiveManeuver/Tstore.mat','Tstore');

%Goal criteria (for the cost function)
% x_g = 0;
% y_g = 10*sin(2*pi*2*Tstore);
theta_g = pi/4;
% 
% xdot_g = 0;
% ydot_g = 10*4*pi*cos(2*pi*2*Tstore); 
% thetadot_g = 0;

%Weighting coefficients
%c1 = xdot, c2 = ydot, c3 = thetadot, c4 = x, c5 = y, c6 = theta
c1 = 1*10^-5; c2 = 1*10^-5; c3 = 10^6; c4 = 10^7; c5 = 10^8; c6 = 10^10; 

%Our initial conditions in the following order:
%Q = [x, y, theta,   phi,  xdot,    ydot, thetadot, phidot]
% q0_og = [0; 0; pi/8; 1.01*pi; 1*10^-4; 1*10^-4; 0; 0];
q0_og = [0, 0, theta_g, (theta_g + pi), 1*10^-4, 1*10^-4, 0, 0];
q0 = q0_og;

%% Import Callin's data
Callin.set1 = csvread('simulationDataset1_blindedX.csv',1);
Callin.set2 = csvread('simulationDataset2_blindedX.csv',1);
Callin.set3 = csvread('simulationDataset3_blindedX.csv',1);

numOfTrajectories(1,1:3) = [size(Callin.set1,1), size(Callin.set2,1),...
    size(Callin.set3,1)];

% ICs = zeros(numOfTrajectories,8);
% ConVars = zeros(numOfTrajectories,3);

%Prescribe the storage data
xx = zeros(timestep*halfwingStrokes,8);
tt = zeros(timestep*halfwingStrokes,1);
% Qstore = zeros(timestep*halfwingStrokes, 8*numOfTrajectories);
% Winstore = zeros(timestep*halfwingStrokes, 8);
% ValSp3 = zeros(numOfTrajectories,12*halfwingStrokes);
% SortedValues = zeros(1,numOfTrajectories); 
% Cost_index = zeros(1,numOfTrajectories);
% ValSp3_ICsave = zeros(halfwingStrokes*8,1);
% costDiagnostic = zeros(9,numOfTrajectories);
% costDiagnostic_2 = zeros(9,wingStrokes);

%% Counters
%Storage vector ranges listed below
%SRC = Store range columns
SRC_1 = 1; SRC_2 = 2; SRC_3 = 3; SRC_4 = 4; 
SRC_5 = 5; SRC_6 = 6; SRC_7 = 7; SRC_8 = 8;

%% Using ode45 

%Start the full run loop
for nn = 1:FullRuns;

ICs = zeros(numOfTrajectories(nn),8);
ConVars = zeros(numOfTrajectories(nn),3);
    
%Unpack Callin's set
    if nn == 1
        A = Callin.set1;
    end

    if nn == 2
        A = Callin.set2;
    end

    if nn == 3
        A = Callin.set3;
    end

Qstore = zeros(timestep*halfwingStrokes, 8*numOfTrajectories(nn));
ValSp3 = zeros(numOfTrajectories(nn),12*halfwingStrokes);
ICs(1:numOfTrajectories(nn),1:8) = A(:,1:8); 
                    % x, y, theta, phi, xdot, ydot, thetadot, phidot
ConVars(1:numOfTrajectories(nn),1:3) = A(:,9:11); % tau, F, alpha

%Start the iteration of consecutive wing strokes
for i = 1:halfwingStrokes;
    
%The storage start and end variables (for storing data in a matrix in the
%correct spot).
SStart = 1 + ((i-1)*timestep);
SEnd = ((i-1)*timestep) + timestep;
VSStart = 1 + ((i-1)*12);
VSEnd = ((i-1)*12) + 12;
NewICStart = VSStart;
NewICEnd = ((i-1)*12) + 8;

% rng('shuffle'); %This re-seeds the random number generator

for j = 1:numOfTrajectories(nn);

q0 = ICs(j,:);

tau0 = ConVars(j,1);
F = ConVars(j,2);
alpha = ConVars(j,3);

thetaR = q0(4) - q0(3) - pi; %This is the resting configuration of our 
    %torsional spring(s) = Initial abdomen angle - initial head angle - pi
% F = 44300*rand(1); %243; %F is the aerodynamic vector in g*cm/(s^2) 
%     %44300 for hawkmoth
% alpha = 2*pi*rand(1); %alpha is the angle of the aerodynamic vector with 
% %respect to the head-thorax in radians
% % tau0 = 10*10000*(2*(rand(1)-0.5)); %The initial torque applied in g*(cm^2)/(s^2)
% % tau0 = 10000*(2*(rand(1)-0.5)); %The initial torque applied in g*(cm^2)/(s^2)
% tau0 = 0; %The initial torque applied in g*(cm^2)/(s^2)
    
OPTIONS = odeset('RelTol',1e-2,'AbsTol',1e-4); 
tic
%[TOUT,YOUT,TE,YE,IE] = ode45(ODEFUN,TSPAN,Y0,OPTIONS)
[T, Q] = ode45(@myODE_3c, [0, 1/hwbf], q0, OPTIONS);

%To clear the pesky warnings/errors
clc;

%tt = linspace(0,T(end),timestep)';
tt = linspace(0,1/hwbf,timestep)';
xx = interp1(T, Q, tt);

%Store ALL of the xx and tt data below. (This will be massive)
% Tstore(SStart:SEnd,1) = [tt + ((i-1)*T(end))];
Qstore(SStart:SEnd,SRC_1) = xx(:,1); %Store the x data
Qstore(SStart:SEnd,SRC_2) = xx(:,2); %Store the y data
Qstore(SStart:SEnd,SRC_3) = xx(:,3); %Store the theta data
Qstore(SStart:SEnd,SRC_4) = xx(:,4); %Store the phi data
Qstore(SStart:SEnd,SRC_5) = xx(:,5); %Store the x dot data
Qstore(SStart:SEnd,SRC_6) = xx(:,6); %Store the y dot data
Qstore(SStart:SEnd,SRC_7) = xx(:,7); %Store the theta dot data
Qstore(SStart:SEnd,SRC_8) = xx(:,8); %Store the phi dot data

%Speed 3 -- The one we primarily care about

%Cost function of speed 3
% Csp3 = c1*(Qstore(SEnd,SRC_5) - xdot_g)^2 +...
%     c2*(Qstore(SEnd,SRC_6) - ydot_g(timestep*i,1))^2 +...
%     c3*(Qstore(SEnd,SRC_7) - thetadot_g)^2 +...
%     c4*(Qstore(SEnd,SRC_1) - x_g)^2 +...
%     c5*(Qstore(SEnd,SRC_2) - y_g(timestep*i,1))^2 +...
%     c6*(Qstore(SEnd,SRC_3) - theta_g)^2;

Csp3 = toc; %Since cost is meaningless in this scenario
    
%Store the control, state vars, and cost which satisfy the conditions
ValSp3(j,VSStart:VSEnd) = [Qstore(SEnd,SRC_1), Qstore(SEnd,SRC_2),...
    Qstore(SEnd,SRC_3), Qstore(SEnd,SRC_4), Qstore(SEnd,SRC_5),...
    Qstore(SEnd,SRC_6), Qstore(SEnd,SRC_7), Qstore(SEnd,SRC_8),...
    tau0, F, alpha, Csp3];

%Now shift our storage vector ranges
SRC_1 = SRC_1+8; SRC_2 = SRC_2+8; SRC_3 = SRC_3+8; SRC_4 = SRC_4+8;
SRC_5 = SRC_5+8; SRC_6 = SRC_6+8; SRC_7 = SRC_7+8; SRC_8 = SRC_8+8;

%Display the progress of the program
disp('Simulation in progress...')
disp(['Full run ', num2str(nn), ' of ', num2str(FullRuns)]);
disp(['Half wingstroke ', num2str(i), ' of ', num2str(halfwingStrokes)]);
disp([num2str(kk), ' of ',...
    num2str(sum(numOfTrajectories)*halfwingStrokes), ' simulations']);
disp([num2str(kkskip-1), ' loops have been skipped']);
disp(['You are ',... 
    num2str(100*kk/(sum(numOfTrajectories)*halfwingStrokes)),...
    '% of the way there!']);

if kk == sum(numOfTrajectories)*halfwingStrokes*FullRuns
    disp('Done!')
end;

%Overall counter (DO NOT RESET)
kk = kk + 1;

end; %End of trajectory simulation (for a given wing stroke)

% Cost_index = find(ValSp3(VSEnd,:)==min(ValSp3(VSEnd,:)));
% 
% if numel(Cost_index) < 1
%     kkskip = kkskip + 1;
%     break 
% end
% 
% %Identify the winners
% WinnerEndConditions = [ValSp3(NewICStart:NewICEnd,Cost_index)]';
% 
% Win_idx_start = find(Qstore(SEnd,:)==WinnerEndConditions(1));
% Win_idx_end = find(Qstore(SEnd,:)==WinnerEndConditions(8));
% 
% if Win_idx_end - Win_idx_start + 1 ~= 8
%     error('Winner indices are not correct')  
% end
% 
% Winstore(SStart:SEnd,1:8) = Qstore(SStart:SEnd, Win_idx_start:Win_idx_end);
% 
% %Sort the values of ValSp3 in order of lowest cost
% % ValSp3_sorted = ValSp3(VSStart:VSEnd,Cost_index);
% ValSp3_sorted = ValSp3(:,Cost_index);
% 
% %To save the ICs of each wing stroke.
% ValSp3_ICsave(VSStart:VSEnd,1)=ValSp3_sorted(VSStart:VSEnd,1);
% 
% %Our NEW initial conditions in the following order:
% %Q = [x, y, theta,   phi,  xdot,    ydot, thetadot, phidot]
% q0(1:8,1) = [ValSp3_sorted(NewICStart:NewICEnd,1)];

%In the event that q0 ends up being all NaNs
if isnan(q0) == 1
    kkskip = kkskip + 1;
    break 
end

%Reset the column counters for the next wing stroke
SRC_1 = 1; SRC_2 = 2; SRC_3 = 3; SRC_4 = 4; 
SRC_5 = 5; SRC_6 = 6; SRC_7 = 7; SRC_8 = 8;

end; %End of entire loop (ALL wing strokes for the full run)

%SAVE the vars necessary for big data stuff
%Four suffixes relevant for our tests:
%am = aggressive maneauver
%con = control
%ntc = control, no tau
%nts = no tau, stiff
%ntf = no tau, floppy
%yts = yes tau, stiff
%ytf = yes tau, floppy

save(['OutputForCallin/Qstore_',num2str(nn+shift),'_am_con.mat'],'Qstore');
% save(['OutputForCallin/Winstore_',num2str(nn+shift),'_am_con.mat'],'Winstore');
save(['OutputForCallin/ValSp3_',num2str(nn+shift),'_am_con.mat'],'ValSp3');
% save(['OutputForCallin/ValSp3_sorted_',num2str(nn+shift),'_am_con.mat'],'ValSp3_sorted');
% save(['OutputForCallin/ValSp3_ICsave_',num2str(nn+shift),'_am_con.mat'],'ValSp3_ICsave');

if kk < numOfTrajectories(nn)*halfwingStrokes*FullRuns && FullRuns > 1
%Now we reset for the next full run
%Our initial conditions in the following order:
%Q = [x, y, theta,   phi,  xdot,    ydot, thetadot, phidot]
q0 = q0_og;

%Prescribe the storage data
xx = zeros(timestep*halfwingStrokes,8);
tt = zeros(timestep*halfwingStrokes,1);
% Tstore = zeros(timestep*wingStrokes, 1);
% Qstore = zeros(timestep*halfwingStrokes, 8*numOfTrajectories);
% ValSp3 = zeros(12*halfwingStrokes, numOfTrajectories);
% SortedValues = zeros(1,numOfTrajectories); 
% Cost_index = zeros(1,numOfTrajectories);
% ValSp3_ICsave = zeros(halfwingStrokes*8,1);
end; 

figure(nn);
hist(ValSp3(:,12)); xlabel('time (seconds)'); ylabel('NumOfSims')

end; %Finished will ALL runs.

disp('Figures WILL NOT render.')
% 
% %% Figures of state variables, trajectories, and flexion w.r.t. time.
% % 
% disp('Figures started rendering.')
% 
% %% Trouble shoot theta
% 
% figure(9);
% for i = 1:halfwingStrokes
%     FigTStart = 1 + ((i-1)*timestep);
%     FigTEnd = ((i-1)*timestep) + timestep;
% plot(Tstore(FigTStart:FigTEnd,1), 180/(pi)*Qstore(FigTStart:FigTEnd,...
%     3:8:end), 'Color', [0.78 0.78 0.78], 'LineWidth',1);
% xlabel('Time (seconds)')
% %ylabel('theta (rad)');
% ylabel('theta (degrees)');
% hold on;
% end
% axis([0, 2, 0, 280])
% 
% %Control variables (Winner only)
% figure(10);
% for i = 1:halfwingStrokes
%     FigTStart = 1 + ((i-1)*timestep);
%     FigTEnd = ((i-1)*timestep) + timestep;
% subplot(5,2,1)
% plot(Tstore(FigTStart:FigTEnd,1), Winstore(FigTStart:FigTEnd,...
%     1), 'r-', 'LineWidth',1);
% % xlabel('Time (seconds)')
% ylabel('x (cm)');
% hold on;
% 
% subplot(5,2,2)
% plot(Tstore(FigTStart:FigTEnd,1), Winstore(FigTStart:FigTEnd,...
%     5), 'r-', 'LineWidth',1);
% % xlabel('Time (seconds)');
% ylabel('x dot (cm/s)');
% hold on;
% 
% subplot(5,2,3)
% plot(Tstore(FigTStart:FigTEnd,1), Winstore(FigTStart:FigTEnd,...
%     2), 'r-', 'LineWidth',1);
% % xlabel('Time (seconds)');
% ylabel('y (cm)');
% hold on;
% 
% subplot(5,2,4)
% plot(Tstore(FigTStart:FigTEnd,1), Winstore(FigTStart:FigTEnd,...
%     6), 'r-', 'LineWidth',1);
% % xlabel('Time (seconds)');
% ylabel('y dot (cm/s)');
% hold on;
% 
% subplot(5,2,5)
% plot(Tstore(FigTStart:FigTEnd,1), (180/pi)*Winstore(FigTStart:FigTEnd,...
%     3), 'r-', 'LineWidth',1);
% % xlabel('Time (seconds)')
% %ylabel('theta (rad)');
% ylabel('theta (degrees)');
% hold on;
% 
% subplot(5,2,6)
% plot(Tstore(FigTStart:FigTEnd,1), (180/pi)*Winstore(FigTStart:FigTEnd,...
%     7), 'r-', 'LineWidth',1);
% %plot(Tstore(:,1), 180/(pi)*Y(:,6), 'LineWidth',2);
% % xlabel('Time (seconds)');
% ylabel('theta dot (deg/s)');
% %ylabel('theta dot (degrees/s)');
% hold on;
% 
% subplot(5,2,7)
% plot(Tstore(FigTStart:FigTEnd,1), (180/pi)*Winstore(FigTStart:FigTEnd,...
%     4), 'r-', 'LineWidth',1);
% % xlabel('Time (seconds)')
% ylabel('phi (degrees)');
% hold on;
% 
% subplot(5,2,8)
% plot(Tstore(FigTStart:FigTEnd,1), (180/pi)*Winstore(FigTStart:FigTEnd,...
%     8), 'r-', 'LineWidth',1);
% %plot(Tstore(:,1), 180/(pi)*Y(:,6), 'LineWidth',2);
% % xlabel('Time (seconds)');
% ylabel('phi dot (deg/s)');
% hold on;
% 
% subplot(5,2,9)
% %plot(Tstore(:,1), Y(:,7), 'LineWidth',2);
% plot(Tstore(FigTStart:FigTEnd,1), (180/pi)*(Winstore(FigTStart:FigTEnd,...
%     4)-Winstore(FigTStart:FigTEnd, 3)-pi),...
%     'r-', 'LineWidth',1);
% xlabel('Time (seconds)');
% ylabel('beta (degrees)');
% hold on;
% 
% subplot(5,2,10)
% plot(Tstore(FigTStart:FigTEnd,1), (180/pi)*(Winstore(FigTStart:FigTEnd,...
%     8)-Winstore(FigTStart:FigTEnd, 7)-pi),...
%     'r-', 'LineWidth',1);
% xlabel('Time (seconds)');
% ylabel('beta dot (deg/s)');
% hold on;
% 
% end;
% 
% disp('  ')
% disp('Figure 10 done')
% 
% %% Our state variables with respect to time
% figure(1);
% for i = 1:halfwingStrokes
%     FigTStart = 1 + ((i-1)*timestep);
%     FigTEnd = ((i-1)*timestep) + timestep;
% subplot(5,2,1)
% plot(Tstore(FigTStart:FigTEnd,1), Qstore(FigTStart:FigTEnd,...
%     1:8:end), 'Color', [0.78 0.78 0.78], 'LineWidth',1);
% xlabel('Time (seconds)')
% ylabel('x (cm)');
% hold on;
% 
% % subplot(5,2,2)
% % plot(Tstore(FigTStart:FigTEnd,1), Qstore(FigTStart:FigTEnd,...
% %     5:8:end), 'Color', [0.78 0.78 0.78], 'LineWidth',1);
% % xlabel('Time (seconds)');
% % ylabel('x dot (cm/s)');
% % hold on;
% 
% subplot(5,2,3)
% plot(Tstore(FigTStart:FigTEnd,1), Qstore(FigTStart:FigTEnd,...
%     2:8:end), 'Color', [0.78 0.78 0.78], 'LineWidth',1);
% xlabel('Time (seconds)');
% ylabel('y (cm)');
% hold on;
% 
% % subplot(5,2,4)
% % plot(Tstore(FigTStart:FigTEnd,1), Qstore(FigTStart:FigTEnd,...
% %     6:8:end), 'Color', [0.78 0.78 0.78], 'LineWidth',1);
% % xlabel('Time (seconds)');
% % ylabel('y dot (cm/s)');
% % hold on;
% 
% subplot(5,2,5)
% %plot(Tstore(:,1), Y(:,5), 'LineWidth',2);
% plot(Tstore(FigTStart:FigTEnd,1), 180/(pi)*Qstore(FigTStart:FigTEnd,...
%     3:8:end), 'Color', [0.78 0.78 0.78], 'LineWidth',1);
% xlabel('Time (seconds)')
% %ylabel('theta (rad)');
% ylabel('theta (degrees)');
% hold on;
% 
% subplot(5,2,6)
% plot(Tstore(FigTStart:FigTEnd,1), 180/(pi)*Qstore(FigTStart:FigTEnd,...
%     7:8:end), 'Color', [0.78 0.78 0.78], 'LineWidth',1);
% %plot(Tstore(:,1), 180/(pi)*Y(:,6), 'LineWidth',2);
% xlabel('Time (seconds)');
% ylabel('theta dot (deg/s)');
% %ylabel('theta dot (degrees/s)');
% hold on;
% 
% subplot(5,2,7)
% %plot(Tstore(:,1), Y(:,5), 'LineWidth',2);
% plot(Tstore(FigTStart:FigTEnd,1), 180/(pi)*Qstore(FigTStart:FigTEnd,...
%     4:8:end), 'Color', [0.78 0.78 0.78], 'LineWidth',1);
% xlabel('Time (seconds)')
% %ylabel('theta (rad)');
% ylabel('phi (degrees)');
% hold on;
% 
% % subplot(5,2,8)
% % plot(Tstore(FigTStart:FigTEnd,1), 180/(pi)*Qstore(FigTStart:FigTEnd,...
% %     8:8:end), 'Color', [0.78 0.78 0.78], 'LineWidth',1);
% % %plot(Tstore(:,1), 180/(pi)*Y(:,6), 'LineWidth',2);
% % xlabel('Time (seconds)');
% % ylabel('phi dot (deg/s)');
% % %ylabel('theta dot (degrees/s)');
% % hold on;
% 
% subplot(5,2,9)
% %plot(Tstore(:,1), Y(:,7), 'LineWidth',2);
% plot(Tstore(FigTStart:FigTEnd,1), 180/(pi)*(Qstore(FigTStart:FigTEnd,...
%     4:8:end)-Qstore(FigTStart:FigTEnd, 3:8:end)-pi),...
%     'Color', [0.78 0.78 0.78], 'LineWidth',1);
% xlabel('Time (seconds)');
% %ylabel('phi (rad)');
% ylabel('beta (degrees)');
% hold on;
% 
% % subplot(5,2,10)
% % %plot(Tstore(:,1), 180/(pi)*Y(:,8), 'LineWidth',2);
% % plot(Tstore(FigTStart:FigTEnd,1), 180/(pi)*(Qstore(FigTStart:FigTEnd,...
% %     8:8:end)-Qstore(FigTStart:FigTEnd, 7:8:end)-pi),...
% %     'Color', [0.78 0.78 0.78], 'LineWidth',1);
% % xlabel('Time (seconds)');
% % ylabel('beta dot (deg/s)');
% % %ylabel('phi dot (degrees/s)');
% % hold on;
% 
% end;
% 
% disp('  ')
% disp('Figure 1 done')
% 
%% All trajectories

AxLmin_x=min(min(Qstore(:,1:8:end))); 
AxLmax_x=max(max(Qstore(:,1:8:end)));
AxLmin_y=min(min(Qstore(:,2:8:end))); 
AxLmax_y=max(max(Qstore(:,2:8:end)));

figure(nn+1);
for i = 1:halfwingStrokes
    FigTStart = 1 + ((i-1)*timestep);
    FigTEnd = ((i-1)*timestep) + timestep;
    plot(Qstore(FigTStart:FigTEnd,1:8:end),...
        Qstore(FigTStart:FigTEnd,2:8:end), 'Color', [0.78 0.78 0.78],...
        'LineWidth',1)
    hold on;
    plot(0,y_g(timestep*i,1),'ko', 'MarkerSize',5)
    xlabel('x (cm)'); ylabel('y (cm)');
    hold on;
end;

for i = 1:halfwingStrokes
    FigTStart = 1 + ((i-1)*timestep);
    FigTEnd = ((i-1)*timestep) + timestep;
    plot(Winstore(FigTStart:FigTEnd,1),...
        Winstore(FigTStart:FigTEnd,2), 'r', 'LineWidth',1)
    hold on;
    plot(Winstore(FigTStart,1), Winstore(FigTStart,2), 'k*')
    hold on;
end;

for i = 1:halfwingStrokes
    plot(0,y_g(timestep*i,1),'ko', 'MarkerSize',5)
    hold on;
end
% disp('Figure 2 done')
% 
% %% All trajectories by wing stroke
% % figure(3);
% % for i = 1:wingStrokes
% %     FigTStart = 1 + ((i-1)*timestep);
% %     FigTEnd = ((i-1)*timestep) + timestep;
% % if wingStrokes > 20 && wingStrokes <= 25
% %     subplot(5,5,i)
% % end
% % if wingStrokes > 15 && wingStrokes <= 20
% %     subplot(4,5,i)
% % end
% % if wingStrokes > 10 && wingStrokes <= 15
% %     subplot(3,5,i) 
% % end
% % if wingStrokes > 5 && wingStrokes <= 10
% %     subplot(2,5,i)
% % end
% % if wingStrokes > 0 && wingStrokes <= 5
% %     subplot(1,wingStrokes,i) 
% % end
% %     plot(Qstore(FigTStart:FigTEnd,1:8:end),...
% %         Qstore(FigTStart:FigTEnd,2:8:end), 'Color', [0.78 0.78 0.78],...
% %         'LineWidth',1);
% %     hold on;
% %     plot(Winstore(FigTStart:FigTEnd,1),...
% %         Winstore(FigTStart:FigTEnd,2), 'r', 'LineWidth',1)
% %     title(['\fontsize{10}Wing stroke ', num2str(i)])
% %     xlabel('x (cm)'); ylabel('y (cm)');
% %     axis([AxLmin_x, AxLmax_x, AxLmin_y, AxLmax_y])
% %     hold on;
% % end;
% % disp('Figure 3 done')
% 
% %% Figures of F, alpha, and tau across all wingstrokes
% 
% %Create a vector that fits ALL of the Force data 
% %Note: These are JUST THE INITIAL CONDITIONS OF CONTROL VARIABLES!
% max_RowLength = size(ValSp3,2);
% Fvec = NaN(max_RowLength,wingStrokes,'double');
% Alphavec = NaN(max_RowLength,wingStrokes,'double');
% Tauvec = NaN(max_RowLength,wingStrokes,'double');
% Costvec = NaN(max_RowLength,wingStrokes,'double');
% 
% %Creating the vector for our state vars
% %Note: These are JUST THE INITIAL CONDITIONS OF STATE VARIABLES!
% xvec = NaN(max_RowLength,wingStrokes,'double');
% yvec = NaN(max_RowLength,wingStrokes,'double');
% thetavec = NaN(max_RowLength,wingStrokes,'double');
% phivec = NaN(max_RowLength,wingStrokes,'double');
% betavec = NaN(max_RowLength,wingStrokes,'double');
% xdotvec = NaN(max_RowLength,wingStrokes,'double');
% ydotvec = NaN(max_RowLength,wingStrokes,'double');
% thetadotvec = NaN(max_RowLength,wingStrokes,'double');
% phidotvec = NaN(max_RowLength,wingStrokes,'double');
% betadotvec = NaN(max_RowLength,wingStrokes,'double');
% 
% for i = 1:wingStrokes
% xvec(1:numel(ValSp3(((i-1)*12+1),:)),i) = ValSp3(((i-1)*12+1),...
%     1:numel(ValSp3(((i-1)*12+1),:)))';
% yvec(1:numel(ValSp3(((i-1)*12+2),:)),i) = ValSp3(((i-1)*12+2),...
%     1:numel(ValSp3(((i-1)*12+2),:)))';
% thetavec(1:numel(ValSp3(((i-1)*12+3),:)),i) = (180/pi)*ValSp3(((i-1)...
%     *12+3), 1:numel(ValSp3(((i-1)*12+3),:)))';
% phivec(1:numel(ValSp3(((i-1)*12+4),:)),i) = (180/pi)*ValSp3(((i-1)...
%     *12+4), 1:numel(ValSp3(((i-1)*12+4),:)))';
% xdotvec(1:numel(ValSp3(((i-1)*12+5),:)),i) = ValSp3(((i-1)*12+5),...
%     1:numel(ValSp3(((i-1)*12+5),:)))';
% ydotvec(1:numel(ValSp3(((i-1)*12+6),:)),i) = ValSp3(((i-1)*12+6),...
%     1:numel(ValSp3(((i-1)*12+6),:)))';
% thetadotvec(1:numel(ValSp3(((i-1)*12+7),:)),i) = (180/pi)*ValSp3(((i-1)...
%     *12+7), 1:numel(ValSp3(((i-1)*12+7),:)))';
% phidotvec(1:numel(ValSp3(((i-1)*12+8),:)),i) = (180/pi)*ValSp3(((i-1)...
%     *12+8), 1:numel(ValSp3(((i-1)*12+8),:)))';
% Fvec(1:numel(ValSp3(((i-1)*12+9),:)),i) = ValSp3(((i-1)*12+9),...
%     1:numel(ValSp3(((i-1)*12+9),:)))';
% Alphavec(1:numel(ValSp3(((i-1)*12+10),:)),i) = (180/pi)*ValSp3(((i-1)*12+10),...
%     1:numel(ValSp3(((i-1)*12+10),:)))';
% Tauvec(1:numel(ValSp3(((i-1)*12+11),:)),i) = ValSp3(((i-1)*12+11),...
%     1:numel(ValSp3(((i-1)*12+11),:)))';
% Costvec(1:numel(ValSp3(((i-1)*12+12),:)),i) = ValSp3(((i-1)*12+12),...
%     1:numel(ValSp3(((i-1)*12+12),:)))';
% end
% 
% %Note the units have already been converted to degrees
% betavec = phivec - thetavec - 180;
% betadotvec = phidotvec - thetadotvec - 180;
% 
% %Fvec = double(Fvec);
% Fvec(Fvec==0) = NaN; Alphavec(Alphavec==0) = NaN; 
% Tauvec(Tauvec==0) = NaN; Costvec(Costvec==0) = NaN;
% xvec(xvec==0) = NaN; yvec(yvec==0) = NaN; 
% thetavec(thetavec==0) = NaN; phivec(phivec==0) = NaN;
% xdotvec(xdotvec==0) = NaN; ydotvec(ydotvec==0) = NaN; 
% thetadotvec(thetadotvec==0) = NaN; phidotvec(phidotvec==0) = NaN;
% 
% figure(4); %Force
% subplot(4,1,1) 
% boxplot([Fvec], 'Notch', 'on')
% ylabel('Force (g*cm*s^-2)');
% hold on;
% 
% %alpha
% subplot(4,1,2)
% boxplot([Alphavec], 'Notch', 'on')
% ylabel('alpha (radians)');
% hold on;
% 
% %tau
% subplot(4,1,3)
% boxplot([Tauvec], 'Notch', 'on')
% ylabel('tau (g*cm^2*s^-2)');
% hold on;
% 
% %Cost
% subplot(4,1,4)
% boxplot([Costvec], 'Notch', 'on')
% ylabel('Cost (unitless)');
% hold on;
% disp('Figure 4 done')
% 
% %% Figures of state variables
% 
% figure(5);
% subplot(5,2,1)
% boxplot([xvec], 'Notch', 'on')
% ylabel('x (cm)');
% hold on;
% 
% subplot(5,2,3)
% boxplot([yvec], 'Notch', 'on')
% ylabel('y (cm)');
% hold on;
% 
% subplot(5,2,5)
% boxplot([thetavec], 'Notch', 'on')
% ylabel('theta (deg)');
% hold on;
% 
% subplot(5,2,7)
% boxplot([phivec], 'Notch', 'on')
% ylabel('phi (deg)');
% hold on;
% 
% subplot(5,2,2)
% boxplot([xdotvec], 'Notch', 'on')
% ylabel('x dot (cm/s)');
% hold on;
% 
% subplot(5,2,4)
% boxplot([ydotvec], 'Notch', 'on')
% ylabel('y dot (cm/s)');
% hold on;
% 
% subplot(5,2,6)
% boxplot([thetadotvec], 'Notch', 'on')
% ylabel('theta dot (deg/s)');
% hold on;
% 
% subplot(5,2,8)
% boxplot([phidotvec], 'Notch', 'on')
% ylabel('phi dot (deg/s)');
% hold on;
% 
% subplot(5,2,9)
% boxplot([betavec], 'Notch', 'on')
% ylabel('beta dot (deg)');
% hold on;
% 
% subplot(5,2,10)
% boxplot([betadotvec], 'Notch', 'on')
% ylabel('beta dot (deg/s)');
% hold on;
% 
% disp('Figure 5 done')
% 
% %% Flexion angle (head-thorax and flexion)
% 
% % Qsp3(Qsp3==0) = NaN;
% % FlexRange = size(Qsp3,2);
% % 
% % figure(#);
% % %Flexion
% % for i = 1:(FlexRange/8);
% % plot(Tstore,((180/(pi))*(Qsp3(:,((i*8)-5)) - Qsp3(:,((i*8)-4))...
% %     - pi)), 'Color', [0.25 0.69 0.61], 'LineWidth',1);
% % xlabel('Time (seconds)');
% % ylabel('Angle (degrees)');
% % legend('Flexion angle');
% % hold on;
% % end;
% % 
% % %Pitch
% % for i = 1:(FlexRange/8);
% % plot(Tstore,((180/(pi))*(Qsp3(:,((i*8)-5)))),...
% %     'Color', [0.17 0.20 0.46],'LineWidth',1, 'DisplayName', 'Pitch angle');
% % xlabel('Time (seconds)');
% % ylabel('Angle (degrees)');
% % hold on;
% % end;
% % 
% % disp('Figure 6 done')
% 
% %% Cost function diagnostics per trajectory WITHOUT total cost 
% 
% % figure(6);
% % plot(1:numOfTrajectories, costDiagnostic(5, :),'.',...
% %     1:numOfTrajectories, costDiagnostic(6, :),'.',...
% %     1:numOfTrajectories, costDiagnostic(7, :),'.',...
% %     1:numOfTrajectories, costDiagnostic(1, :),'x',...
% %     1:numOfTrajectories, costDiagnostic(2, :),'.',...
% %     1:numOfTrajectories, costDiagnostic(3, :),'.');
% % xlabel('Trajectory number');
% % ylabel('Cost (unitless)');
% % legend('x dot','y dot','theta dot','x','y','theta');
% % disp('Figure 6 done')
% % 
% % %% Cost function diagnostics per trajectory WITH total cost 
% % 
% % figure(7);
% % plot(1:numOfTrajectories, costDiagnostic(5, :),'.',...
% %     1:numOfTrajectories, costDiagnostic(6, :),'.',...
% %     1:numOfTrajectories, costDiagnostic(7, :),'.',...
% %     1:numOfTrajectories, costDiagnostic(1, :),'x',...
% %     1:numOfTrajectories, costDiagnostic(2, :),'.',...
% %     1:numOfTrajectories, costDiagnostic(3, :),'.',...
% %     1:numOfTrajectories, costDiagnostic(9,:), 'o', 'MarkerSize',4);
% % xlabel('Trajectory number');
% % ylabel('Cost (unitless)');
% % legend('x dot','y dot','theta dot','x','y','theta','Total cost');
% % disp('Figure 7 done')
% % 
% % %% Cost function diagnostics per trajectory (likely not gonna be useful)
% % 
% % figure(8);
% % plot(1:wingStrokes, costDiagnostic_2(5, :),...
% %     1:wingStrokes, costDiagnostic_2(6, :),...
% %     1:wingStrokes, costDiagnostic_2(7, :),...
% %     1:wingStrokes, costDiagnostic_2(1, :),...
% %     1:wingStrokes, costDiagnostic_2(2, :),...
% %     1:wingStrokes, costDiagnostic_2(3, :), 'LineWidth', 1);
% % xlabel('Wing stroke');
% % ylabel('Cost (unitless)');
% % legend('x dot','y dot','theta dot','x','y','theta');
% % disp('Figure 8 done')

%% Final display on window
disp('   ');
disp('Completely done, son.');