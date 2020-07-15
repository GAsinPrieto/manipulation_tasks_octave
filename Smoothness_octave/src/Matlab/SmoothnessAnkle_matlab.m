%% Smoothness of ankle angular velocity signal: SPARC, NP an LDLJ

% Code AUTHOR: Yaiza Benito Molpeceres. DATE: January-May 2020.

% Proposed only for Sagittal Lifting tasks
% Calculated for two phases named:
    % Lifting phase: Subject holds and deposits the box (idxOLF:idxFLF)
    % Lowering phase: Subject holds and deposits the box in its original location (idxOLW:idxFLW)

clear all % Clear variables
close all % Close figures
clc

load('..\tests\data\input\dinamica57_B.mat')

Ts = 1/Fs;
parameters = [0.025; Fs/2; 4];  % Threshold 0.025 or less, between [0.015, 0.025]

%% Import model data for Ankle Relative Angles

% Ankle: 1 dof: flexo-extension

        RAnkleAngles_x = ModelData.Raw.(ModelOutput{5})(1,:); 
        RAnkleAngles_y = ModelData.Raw.(ModelOutput{5})(2,:);
        RAnkleAngles_z = ModelData.Raw.(ModelOutput{5})(3,:);
        

        LAnkleAngles_x = ModelData.Raw.(ModelOutput{6})(1,:); 
        LAnkleAngles_y = ModelData.Raw.(ModelOutput{6})(2,:);
        LAnkleAngles_z = ModelData.Raw.(ModelOutput{6})(3,:);
         
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
plot(t, position_Z, t, RAnkleAngles_x, t, RAnkleAngles_y, t, RAnkleAngles_z, t, LAnkleAngles_x, t, LAnkleAngles_y, t, LAnkleAngles_z);
legend('Z', 'RAnkleAngles_x', 'RAnkleAngles_y', 'RAnkleAngles_z', 'LAnkleAngles_x', 'LAnkleAngles_y', 'LAnkleAngles_z');
xlabel('Time in seconds');
title('Vertical Box Trajectory and Ankle Relative Angles');


%% Angular Velocities Calculation(rad/s)

    RAnkle_angVelocity_x = dv(t,RAnkleAngles_x);
    RAnkle_angVelocity_x = RAnkle_angVelocity_x'/Ts;
    RAnkle_angVelocity_y = dv(t,RAnkleAngles_y);
    RAnkle_angVelocity_y = RAnkle_angVelocity_y'/Ts;
    RAnkle_angVelocity_z = dv(t,RAnkleAngles_z);
    RAnkle_angVelocity_z = RAnkle_angVelocity_z'/Ts;

    LAnkle_angVelocity_x = dv(t,LAnkleAngles_x);
    LAnkle_angVelocity_x = LAnkle_angVelocity_x'/Ts;
    LAnkle_angVelocity_y = dv(t,LAnkleAngles_y);
    LAnkle_angVelocity_y = LAnkle_angVelocity_y'/Ts;
    LAnkle_angVelocity_z = dv(t,LAnkleAngles_z);
    LAnkle_angVelocity_z = LAnkle_angVelocity_z'/Ts;  
    
%% Crop Relative angular velocities to lifting phase
       
    RAnkle_angVelocity_x_lifting_sparc =  RAnkle_angVelocity_x(idxOLF:idxFLF);
    RAnkle_angVelocity_y_lifting_sparc =  RAnkle_angVelocity_y(idxOLF:idxFLF);
    RAnkle_angVelocity_z_lifting_sparc =  RAnkle_angVelocity_z(idxOLF:idxFLF);

    LAnkle_angVelocity_x_lifting_sparc =  LAnkle_angVelocity_x(idxOLF:idxFLF);
    LAnkle_angVelocity_y_lifting_sparc =  LAnkle_angVelocity_y(idxOLF:idxFLF);
    LAnkle_angVelocity_z_lifting_sparc =  LAnkle_angVelocity_z(idxOLF:idxFLF);
    
figure (3) 
t2 = (0:1:(length(RAnkle_angVelocity_x_lifting_sparc)-1))*Ts;
subplot (2,1,1) 
plot(t2, RAnkle_angVelocity_x_lifting_sparc, t2, RAnkle_angVelocity_y_lifting_sparc, t2, RAnkle_angVelocity_z_lifting_sparc);
legend('x', 'y', 'z');
xlabel('Time in seconds');
ylabel('Angular velocity in degrees/s');
title ('Cropped Right Ankle Angular Velocities Lifting Phase')
subplot (2,1,2)
plot(t2, LAnkle_angVelocity_x_lifting_sparc, t2, LAnkle_angVelocity_y_lifting_sparc, t2, LAnkle_angVelocity_z_lifting_sparc);
legend('x', 'y', 'z');
xlabel('Time in seconds');
ylabel('Angular velocity in degrees/s');
title ('Cropped Left Ankle Angular Velocities Lifting Phase')


%% Smoothness Lifting Ankle: SPARC, NP and LDLJ

S_lifting_Ankle = {'S_R_lifting_X', 'S_R_lifting_Y', 'S_R_lifting_Z', 'S_L_lifting_X', 'S_L_lifting_Y', 'S_L_lifting_Z'};
    Ankle_angVel(1,:) = RAnkle_angVelocity_x_lifting_sparc;
    Ankle_angVel(2,:) = RAnkle_angVelocity_y_lifting_sparc;
    Ankle_angVel(3,:) = RAnkle_angVelocity_z_lifting_sparc;
    Ankle_angVel(4,:) = LAnkle_angVelocity_x_lifting_sparc;
    Ankle_angVel(5,:) = LAnkle_angVelocity_y_lifting_sparc;
    Ankle_angVel(6,:) = LAnkle_angVelocity_z_lifting_sparc;

for i = 1:length(S_lifting_Ankle)
    
    speed = Ankle_angVel(i,:);
    %SPARC
    S_sparc = SpectralArcLength(speed, Ts, parameters);
    [S_lifting_Ankle_SPARC{i}] = S_sparc;
    
    %NP
    [S_np, locs_np] = np(t2,speed);
    [S_lifting_Ankle_NP{i}] = S_np;
    [peaks_lifting{i}] = locs_np;
    
    %LDLJ
    [Dlj S_ldlj] = ldlj(t2,speed);
    [S_lifting_Ankle_LDLJ{i}] = S_ldlj;
end 



%% Events of interest in vertical box trajectory identification

% Lowering phase

flat3 = idx < locs(2);
idx_flat3 = find(flat3);

flat4 = idx > locs(2);
idx_flat4 = find(flat4);

idxOLW = idx(idx_flat3(length(idx_flat3)));
idxFLW = idx(idx_flat4(3));                                       % para din 59 es 4, para din 03 y 61 es 2 y para el resto es 1

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
      
    
    RAnkle_angVelocity_x_lowering_sparc =  RAnkle_angVelocity_x(idxOLW:idxFLW);
    RAnkle_angVelocity_y_lowering_sparc =  RAnkle_angVelocity_y(idxOLW:idxFLW);
    RAnkle_angVelocity_z_lowering_sparc =  RAnkle_angVelocity_z(idxOLW:idxFLW);

    LAnkle_angVelocity_x_lowering_sparc =  LAnkle_angVelocity_x(idxOLW:idxFLW);
    LAnkle_angVelocity_y_lowering_sparc =  LAnkle_angVelocity_y(idxOLW:idxFLW);
    LAnkle_angVelocity_z_lowering_sparc =  LAnkle_angVelocity_z(idxOLW:idxFLW);
    
figure (5) 
t3 = (0:1:(length(RAnkle_angVelocity_x_lowering_sparc)-1))*Ts;
subplot (2,1,1)
plot(t3, RAnkle_angVelocity_x_lowering_sparc, t3, RAnkle_angVelocity_y_lowering_sparc, t3, RAnkle_angVelocity_z_lowering_sparc);
legend('x', 'y', 'z');
xlabel('Time in seconds');
ylabel('Angular velocity in degrees/s');
title ('Cropped Right Ankle Angular Velocities Lifting Phase')
subplot (2,1,2)
plot(t3, LAnkle_angVelocity_x_lowering_sparc, t3, LAnkle_angVelocity_y_lowering_sparc, t3, LAnkle_angVelocity_z_lowering_sparc);
legend('x', 'y', 'z');
xlabel('Time in seconds');
ylabel('Angular velocity in degrees/s');
title ('Cropped Left Ankle Angular Velocities Lowering Phase')

%% Smoothness Lowering Ankle: SPARC, NP and LDLJ

S_lowering_Ankle = {'S_R_lowering_X', 'S_R_lowering_Y', 'S_R_lowering_Z', 'S_L_lowering_X', 'S_L_lowering_Y', 'S_L_lowering_Z'};
    Ankle_angVel_low(1,:) =  RAnkle_angVelocity_x_lowering_sparc;
    Ankle_angVel_low(2,:) =  RAnkle_angVelocity_y_lowering_sparc;
    Ankle_angVel_low(3,:) =  RAnkle_angVelocity_z_lowering_sparc;
    Ankle_angVel_low(4,:) =  LAnkle_angVelocity_x_lowering_sparc;
    Ankle_angVel_low(5,:) =  LAnkle_angVelocity_y_lowering_sparc;
    Ankle_angVel_low(6,:) =  LAnkle_angVelocity_z_lowering_sparc;

for i = 1:length(S_lowering_Ankle)
    
    speed = Ankle_angVel_low(i,:);
    %SPARC
     S_sparc = SpectralArcLength(speed, Ts, parameters);
    [S_lowering_Ankle_SPARC{i}] = S_sparc;
    
    %NP
    [S_np, locs_np] = np(t3,speed);
    [S_lowering_Ankle_NP{i}] = S_np;
    [peaks_lowering{i}] = locs_np;
    
    %LDLJ
    [Dlj S_ldlj] = ldlj(t3,speed);
    [S_lowering_Ankle_LDLJ{i}] = S_ldlj;
end 
   

Smoothness_Lifting_Ankle = table([S_lifting_Ankle_SPARC(1); S_lifting_Ankle_SPARC(4); S_lifting_Ankle_SPARC(2); S_lifting_Ankle_SPARC(5);S_lifting_Ankle_SPARC(3);S_lifting_Ankle_SPARC(6)],[S_lifting_Ankle_NP(1); S_lifting_Ankle_NP(4); S_lifting_Ankle_NP(2); S_lifting_Ankle_NP(5);S_lifting_Ankle_NP(3);S_lifting_Ankle_NP(6)],[S_lifting_Ankle_LDLJ(1); S_lifting_Ankle_LDLJ(4); S_lifting_Ankle_LDLJ(2); S_lifting_Ankle_LDLJ(5);S_lifting_Ankle_LDLJ(3);S_lifting_Ankle_LDLJ(6)],'VariableNames',{'SPARC','NP','LDLJ'},'RowNames',{'Right X','Left  X','Right Y','Left  Y','Right Z','Left  Z'})
Smoothness_Lowering_Ankle = table([S_lowering_Ankle_SPARC(1); S_lowering_Ankle_SPARC(4); S_lowering_Ankle_SPARC(2); S_lowering_Ankle_SPARC(5);S_lowering_Ankle_SPARC(3);S_lowering_Ankle_SPARC(6)],[S_lowering_Ankle_NP(1); S_lowering_Ankle_NP(4); S_lowering_Ankle_NP(2); S_lowering_Ankle_NP(5);S_lowering_Ankle_NP(3);S_lowering_Ankle_NP(6)],[S_lowering_Ankle_LDLJ(1); S_lowering_Ankle_LDLJ(4); S_lowering_Ankle_LDLJ(2); S_lowering_Ankle_LDLJ(5);S_lowering_Ankle_LDLJ(3);S_lowering_Ankle_LDLJ(6)],'VariableNames',{'SPARC','NP','LDLJ'},'RowNames',{'Right X','Left  X','Right Y','Left  Y','Right Z','Left  Z'})

