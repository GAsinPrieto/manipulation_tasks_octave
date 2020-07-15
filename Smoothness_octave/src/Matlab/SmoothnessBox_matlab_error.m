%% Smoothness of box trajectory marker

% Code AUTHOR: Yaiza Benito Molpeceres. DATE: January-May 2020.

% Proposed only for Sagittal Lifting tasks
% It is observed in Vicon captures that when the subject deposits a loaded box 
% it bounces. Thus, measuring its smoothness can show 
% differences regarding precision or delicacy between placing a loaded box
% and a unloaded one. 

clear all % Clear variables
close all % Close figures
clc

load('..\tests\data\input\dinamica59_B.mat')

Ts = 1/Fs;
parameters = [0.015; Fs/2; 4];   % threshold between [0.015,0.020]


%% Vertical SIDEBOX1 box marker trajectory

    position_Z = SIDEBOX1(:,3)';
    t = (0:1:(length(position_Z)-1))*Ts;
    
%% Angular Velocity Calculation(rad/s)

    Box_angVelocity_z = dv(t,position_Z);
    Box_angVelocity_z = Box_angVelocity_z'/Ts;
    
    
%% Events of interest in vertical box trajectory identification
 
position_Z = SIDEBOX3(:,3)';
 
[pks, locs] = findpeaks(position_Z,'minPeakProminence',4);
 
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

% Lowering phase

flat3 = idx < locs(2);
idx_flat3 = find(flat3);

flat4 = idx > locs(2);
idx_flat4 = find(flat4);

idxOLW = idx(idx_flat3(length(idx_flat3)));
idxFLW = idx(idx_flat4(3));

figure(2)
plot(t,position_Z,t(locs(2)),position_Z(locs(2)),'o');
hold on
plot(t,position_Z,'r*','MarkerIndices',idxOLW );
hold on
plot(t,position_Z,'r*','MarkerIndices',idxFLW);
hold off
xlabel('Time in seconds');
ylabel('Millimetres');
title('Events of interest in box trajectory during lowering for SPARC');


%% Crop relative angular velocities to lifting phase

    Box_angVelocity_z_lifting_sparc =  Box_angVelocity_z(idxOLF:idxFLF);
    
    Box_angVelocity_z_lowering_sparc =  Box_angVelocity_z(idxOLW:idxFLW);

    %t2 = (0:1:(length(Box_angVelocity_z_lifting_sparc)-1))*Ts;
    
    %t3 = (0:1:(length(Box_angVelocity_z_lowering_sparc)-1))*Ts;
    
%% SPARC algorithm: smoothness calculation

    speed_lifting = Box_angVelocity_z_lifting_sparc(1,:);
    S_lifting_SPARC = SpectralArcLength(speed_lifting, Ts, parameters);
    
    speed_lowering = Box_angVelocity_z_lowering_sparc(1,:);
    S_lowering_SPARC = SpectralArcLength(speed_lowering, Ts, parameters);
    
%% NP
    
    speed_lf_np = Box_angVelocity_z_lifting_sparc(1,:);
    [S_lifting_NP, locs_np] = np(t2,speed_lf_np);
   
    speed_lw_np = Box_angVelocity_z_lowering_sparc(1,:);
    [S_lowering_NP, locs2_np] = np(t3,speed_lw_np);

%% LDLJ
    
    speed_lf_ldlj = Box_angVelocity_z_lifting_sparc(1,:);
    [Dlj S_lifting_LDLJ] = ldlj(t2,speed_lf_ldlj);
       
    speed_lw_ldlj = Box_angVelocity_z_lowering_sparc(1,:);
    [Dlj S_lowering_LDLJ] = ldlj(t3,speed_lw_ldlj);


%% Table

Smoothness_Lifting_Box = table([S_lifting_SPARC],[S_lifting_NP],[S_lifting_LDLJ],'VariableNames',{'SPARC','NP','LDLJ'},'RowNames',{'Box Z'})
Smoothness_Lowering_Box = table([S_lowering_SPARC],[S_lowering_NP],[S_lowering_LDLJ],'VariableNames',{'SPARC','NP','LDLJ'},'RowNames',{'Box Z'})

