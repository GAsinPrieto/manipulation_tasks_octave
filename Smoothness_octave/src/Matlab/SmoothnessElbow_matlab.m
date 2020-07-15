%% Smoothness of elbow angular velocity signal: SPARC, NP an LDLJ

% Code AUTHOR: Yaiza Benito Molpeceres. DATE: January-May 2020.

% Proposed only for Sagittal Lifting tasks
% Calculated for two phases named:
    % Lifting phase: Subject holds and deposits the box (idxOLF:idxFLF)
    % Lowering phase: Subject holds and deposits the box in its original location (idxOLW:idxFLW)

clear all % Clear variables
close all % Close figures
clc

load('../tests/data/input/dinamica59_B.mat')

Ts = 1/Fs;
parameters = [0.025; Fs/2; 4];    % Threshold 0.025 or less, between [0.015, 0.025]

%% Import model data for Elbow Relative Angles

% Elbow: 1 dof flexion-extension

        RElbowAngles_x = ModelData.Raw.(ModelOutput{7})(1,:); 
        RElbowAngles_y = ModelData.Raw.(ModelOutput{7})(2,:);
        RElbowAngles_z = ModelData.Raw.(ModelOutput{7})(3,:);
        

        LElbowAngles_x = ModelData.Raw.(ModelOutput{8})(1,:);
        LElbowAngles_y = ModelData.Raw.(ModelOutput{8})(2,:);
        LElbowAngles_z = ModelData.Raw.(ModelOutput{8})(3,:);
        
%% Events of interest in vertical box trajectory identification
 
position_Z = SIDEBOX3(:,3)';
 
[pks, locs] = findpeaks(position_Z,'minPeakProminence',4);
t = (0:1:(length(position_Z)-1))*Ts;
 
TF1 = islocalmin(position_Z, 'FlatSelection','first');
%for i = 1:frames                               % solo para dinamica 03
 %   TF2(1,i) = position_Z(1,i) > 0;
%end
idx = find(TF1);
%idx_03 = find(TF2);
flat = idx < locs(1);
idx_flat = find(flat);
flat2 = idx > locs(1);
idx_flat2 = find(flat2);
 
flat3 = idx < locs(2);
idx_flat3 = find(flat3);
 
flat4 = idx > locs(2);
idx_flat4 = find(flat4);
 
% Indexes definiton for signal segmentation
idxO1 = 1;
idxF1 = idx(idx_flat(length(idx_flat)));
 
 
% Indexes definiton for signal segmentation
idxOLF = idx(idx_flat(length(idx_flat)));
idxFLF = idx(idx_flat2(1));


figure(1)
plot(t,position_Z,t(locs(1)),position_Z(locs(1)),'o');
hold on
plot(t,position_Z,'r*','MarkerIndices',idxOLF);
hold on
plot(t,position_Z,'r*','MarkerIndices',idxFLF);
hold off
xlabel('Time in seconds');
ylabel('Millimetres');
title('Events of interest in box trajectory during lifting for SPARC');

figure(2)
plot(t, position_Z, t, RElbowAngles_x, t, RElbowAngles_y, t, RElbowAngles_z, t, LElbowAngles_x, t, LElbowAngles_y, t, LElbowAngles_z);
legend('Z', 'RElbowAngles_x', 'RElbowAngles_y', 'RElbowAngles_z', 'LElbowAngles_x', 'LElbowAngles_y', 'LElbowAngles_z');
xlabel('Time in seconds');
title('Vertical Box Trajectory and Elbow Relative Angles');


%% Angular Velocities Calculation(rad/s)

    RElbow_angVelocity_x = dv(t,RElbowAngles_x);
    RElbow_angVelocity_x = RElbow_angVelocity_x'/Ts;
    RElbow_angVelocity_y = dv(t,RElbowAngles_y);
    RElbow_angVelocity_y = RElbow_angVelocity_y'/Ts;
    RElbow_angVelocity_z = dv(t,RElbowAngles_z);
    RElbow_angVelocity_z = RElbow_angVelocity_z'/Ts;

    LElbow_angVelocity_x = dv(t,LElbowAngles_x);
    LElbow_angVelocity_x = LElbow_angVelocity_x'/Ts;
    LElbow_angVelocity_y = dv(t,LElbowAngles_y);
    LElbow_angVelocity_y = LElbow_angVelocity_y'/Ts;
    LElbow_angVelocity_z = dv(t,LElbowAngles_z);
    LElbow_angVelocity_z = LElbow_angVelocity_z'/Ts;  
      

%% Crop Relative angular velocities to lifting phase

    RElbow_angVelocity_x_lifting_sparc =  RElbow_angVelocity_x(idxOLF:idxFLF);
    RElbow_angVelocity_y_lifting_sparc =  RElbow_angVelocity_y(idxOLF:idxFLF);
    RElbow_angVelocity_z_lifting_sparc =  RElbow_angVelocity_z(idxOLF:idxFLF);

    LElbow_angVelocity_x_lifting_sparc =  LElbow_angVelocity_x(idxOLF:idxFLF);
    LElbow_angVelocity_y_lifting_sparc =  LElbow_angVelocity_y(idxOLF:idxFLF);
    LElbow_angVelocity_z_lifting_sparc =  LElbow_angVelocity_z(idxOLF:idxFLF);
    
figure (3) 
t2 = (0:1:(length(RElbow_angVelocity_x_lifting_sparc)-1))*Ts;
subplot (2,1,1)
plot(t2, RElbow_angVelocity_x_lifting_sparc, t2, RElbow_angVelocity_y_lifting_sparc, t2, RElbow_angVelocity_z_lifting_sparc);
legend('x', 'y', 'z');
xlabel('Time in seconds');
ylabel('Angular velocity in degrees/s');
title ('Cropped Right Elbow Angular Velocities Lifting Phase')
subplot (2,1,2)
plot(t2, LElbow_angVelocity_x_lifting_sparc, t2, LElbow_angVelocity_y_lifting_sparc, t2, LElbow_angVelocity_z_lifting_sparc);
legend('x', 'y', 'z');
xlabel('Time in seconds');
ylabel('Angular velocity in degrees/s');
title ('Cropped Left Elbow Angular Velocities Lifting Phase')

%% Smoothness Lifting Elbow: SPARC, NP and LDLJ

S_lifting_Elbow = {'S_R_lifting_X', 'S_R_lifting_Y', 'S_R_lifting_Z', 'S_L_lifting_X', 'S_L_lifting_Y', 'S_L_lifting_Z'};
    Elbow_angVel(1,:) = RElbow_angVelocity_x_lifting_sparc;
    Elbow_angVel(2,:) = RElbow_angVelocity_y_lifting_sparc;
    Elbow_angVel(3,:) = RElbow_angVelocity_z_lifting_sparc;
    Elbow_angVel(4,:) = LElbow_angVelocity_x_lifting_sparc;
    Elbow_angVel(5,:) = LElbow_angVelocity_y_lifting_sparc;
    Elbow_angVel(6,:) = LElbow_angVelocity_z_lifting_sparc;

for i = 1:length(S_lifting_Elbow)
    
    speed = Elbow_angVel(i,:);
    %SPARC
    S_sparc = SpectralArcLength(speed, Ts, parameters);
    [S_lifting_Elbow_SPARC{i}] = S_sparc;
    
    %NP
    [S_np, locs_np] = np(t2,speed);
    [S_lifting_Elbow_NP{i}] = S_np;
    [peaks_lifting{i}] = locs_np;
    
    %LDLJ
    [Dlj S_ldlj] = ldlj(t2,speed);
    [S_lifting_Elbow_LDLJ{i}] = S_ldlj;
end 


%% Events of interest in vertical box trajectory identification

% Lowering phase

flat3 = idx < locs(2);
idx_flat3 = find(flat3);

flat4 = idx > locs(2);
idx_flat4 = find(flat4);

idxOLW = idx(idx_flat3(length(idx_flat3)));
idxFLW = idx(idx_flat4(3));                                 % para din 59 es 4, para din 03 y 61 es 2 y para el resto es 1

figure(4)
plot(t,position_Z,t(locs(2)),position_Z(locs(2)),'o');
hold on
plot(t,position_Z,'r*','MarkerIndices',idxOLW);
hold on
plot(t,position_Z,'r*','MarkerIndices',idxFLW);
hold off
xlabel('Time in seconds');
ylabel('Millimetres');
title('Events of interest in box trajectory during lowering for SPARC');


%% Crop Relative angular velocities to lowering phase


    RElbow_angVelocity_x_lowering_sparc =  RElbow_angVelocity_x(idxOLW:idxFLW);
    RElbow_angVelocity_y_lowering_sparc =  RElbow_angVelocity_y(idxOLW:idxFLW);
    RElbow_angVelocity_z_lowering_sparc =  RElbow_angVelocity_z(idxOLW:idxFLW);

    LElbow_angVelocity_x_lowering_sparc =  LElbow_angVelocity_x(idxOLW:idxFLW);
    LElbow_angVelocity_y_lowering_sparc =  LElbow_angVelocity_y(idxOLW:idxFLW);
    LElbow_angVelocity_z_lowering_sparc =  LElbow_angVelocity_z(idxOLW:idxFLW);
    
figure (5) 
t3 = (0:1:(length(RElbow_angVelocity_x_lowering_sparc)-1))*Ts;
subplot (2,1,1)
plot(t3, RElbow_angVelocity_x_lowering_sparc, t3, RElbow_angVelocity_y_lowering_sparc, t3, RElbow_angVelocity_z_lowering_sparc);
legend('x', 'y', 'z');
xlabel('Time in seconds');
ylabel('Angular velocity in degrees/s');
title ('Cropped Right Elbow Angular Velocities Lowering Phase')
subplot (2,1,2)
plot(t3, LElbow_angVelocity_x_lowering_sparc, t3, LElbow_angVelocity_y_lowering_sparc, t3, LElbow_angVelocity_z_lowering_sparc);
legend('x', 'y', 'z');
xlabel('Time in seconds');
ylabel('Angular velocity in degrees/s');
title ('Cropped Left Elbow Angular Velocities Lowering Phase')

%% Smoothness Lowering Elbow: SPARC, NP and LDLJ

S_lowering_Elbow = {'S_R_lowering_X', 'S_R_lowering_Y', 'S_R_lowering_Z', 'S_L_lowering_X', 'S_L_lowering_Y', 'S_L_lowering_Z'};
    Elbow_angVel_low(1,:) =  RElbow_angVelocity_x_lowering_sparc;
    Elbow_angVel_low(2,:) =  RElbow_angVelocity_y_lowering_sparc;
    Elbow_angVel_low(3,:) =  RElbow_angVelocity_z_lowering_sparc;
    Elbow_angVel_low(4,:) =  LElbow_angVelocity_x_lowering_sparc;
    Elbow_angVel_low(5,:) =  LElbow_angVelocity_y_lowering_sparc;
    Elbow_angVel_low(6,:) =  LElbow_angVelocity_z_lowering_sparc;

for i = 1:length(S_lowering_Elbow)
    
    speed = Elbow_angVel_low(i,:);
    %SPARC
     S_sparc = SpectralArcLength(speed, Ts, parameters);
    [S_lowering_Elbow_SPARC{i}] = S_sparc;
    
    %NP
    [S_np, locs_np] = np(t3,speed);
    [S_lowering_Elbow_NP{i}] = S_np;
    [peaks_lowering{i}] = locs_np;
    
    %LDLJ
    [Dlj S_ldlj] = ldlj(t3,speed);
    [S_lowering_Elbow_LDLJ{i}] = S_ldlj;
end 
  

Smoothness_Lifting_Elbow = table([S_lifting_Elbow_SPARC(1); S_lifting_Elbow_SPARC(4); S_lifting_Elbow_SPARC(2); S_lifting_Elbow_SPARC(5);S_lifting_Elbow_SPARC(3);S_lifting_Elbow_SPARC(6)],[S_lifting_Elbow_NP(1); S_lifting_Elbow_NP(4); S_lifting_Elbow_NP(2); S_lifting_Elbow_NP(5);S_lifting_Elbow_NP(3);S_lifting_Elbow_NP(6)],[S_lifting_Elbow_LDLJ(1); S_lifting_Elbow_LDLJ(4); S_lifting_Elbow_LDLJ(2); S_lifting_Elbow_LDLJ(5);S_lifting_Elbow_LDLJ(3);S_lifting_Elbow_LDLJ(6)],'VariableNames',{'SPARC','NP','LDLJ'},'RowNames',{'Right X','Left  X','Right Y','Left  Y','Right Z','Left  Z'})
Smoothness_Lowering_Elbow = table([S_lowering_Elbow_SPARC(1); S_lowering_Elbow_SPARC(4); S_lowering_Elbow_SPARC(2); S_lowering_Elbow_SPARC(5);S_lowering_Elbow_SPARC(3);S_lowering_Elbow_SPARC(6)],[S_lowering_Elbow_NP(1); S_lowering_Elbow_NP(4); S_lowering_Elbow_NP(2); S_lowering_Elbow_NP(5);S_lowering_Elbow_NP(3);S_lowering_Elbow_NP(6)],[S_lowering_Elbow_LDLJ(1); S_lowering_Elbow_LDLJ(4); S_lowering_Elbow_LDLJ(2); S_lowering_Elbow_LDLJ(5);S_lowering_Elbow_LDLJ(3);S_lowering_Elbow_LDLJ(6)],'VariableNames',{'SPARC','NP','LDLJ'},'RowNames',{'Right X','Left  X','Right Y','Left  Y','Right Z','Left  Z'})

