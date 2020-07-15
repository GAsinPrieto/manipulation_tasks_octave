%% Total time to perform the task in Lateral Box Transfer (exoskeleton segmentation)

% Code AUTHOR: Yaiza Benito Molpeceres. DATE: January-May 2020.

clear all % Clear variables
close all % Close figures
clc

load('../tests/data/input/dinamica56_B.mat')
             
Ts = 1/Fs;

%% Total time to perform the task 

t_total = (double(frames)*Ts);

%% Events of interest in feet markers trajectories

% Feet markers for signal segmentation
LTOE_z = LTOE(:,3)'; 
RTOE_z = RTOE(:,3)';   

% Left toe:
 
[pks, locs] = findpeaks(LTOE_z,'minPeakProminence',10,'MinPeakHeight',80);                    
   
    
TF1 = islocalmin(LTOE_z, 'FlatSelection', 'first');

idx = find(TF1);
flat = idx < locs(1);
idx_flat = find(flat);
 
flat2 = idx < locs(3)& idx > locs(2);
idx_flat2 = find(flat2);
 
flat3 = idx > locs(4);
idx_flat3 = find(flat3);
 
 
% Right Toe: 
 
[pks2, locs2] = findpeaks(RTOE_z,'minPeakProminence',10,'MinPeakHeight',80);  
 
TF2 = islocalmin(RTOE_z, 'FlatSelection', 'first');

idx2 = find(TF2);
flat4 = idx2 < locs2(1);
idx_flat4 = find(flat4);
 
flat5 = idx2 > locs2(length(locs2));
idx_flat5 = find(flat5);
flat6 = idx2 > locs2(3);
idx_flat6 = find(flat6);
 
 
% Indexes definiton for signal segmentation
idxO1 = 1;
idxF1 = idx2(idx_flat4(length(idx_flat4))); %idxO2
idxF2 = idx2(idx_flat6(1)); %idxO3
idxF3 = idx(idx_flat2(length(idx_flat2))); %idxO4
idxF4 = idx(idx_flat3(4)); %idxO5
idxF5 = length(LTOE_z);
 
figure(1)
t= 0:(length(LTOE)-1);
plot(t, zscore(LTOE_z));
hold on
plot(t,zscore(RTOE_z));
plot(t,zscore(LTOE_z),'r*', 'MarkerIndices', idxO1);
plot(t,zscore(RTOE_z),'r*', 'MarkerIndices', idxF2);
plot(t,zscore(LTOE_z),'r*', 'MarkerIndices', idxF3);
plot(t,zscore(LTOE_z),'r*', 'MarkerIndices', idxF5);
plot(t,zscore(RTOE_z),'r*', 'MarkerIndices', idxO1);
plot(t,zscore(RTOE_z),'r*', 'MarkerIndices', idxF1);
plot(t,zscore(LTOE_z),'r*', 'MarkerIndices', idxF4);
plot(t,zscore(RTOE_z),'r*', 'MarkerIndices', idxF5);
hold off
legend('Left Toe Z', 'Right Toe Z');
title('Segmentation Exoskeleton Trials in 5 phases')


% Five phases shown in figure 1:
    % 1: Subject picks up the box in the sagittal plane (idxO1:idxF1)
    % 2: Subject takes some steps to rotate to the frontal plane (idxF1:idxF2)
    % 3: Subject deposits the box in the frontal plane (idxF2:idxF3)
    % 4: Subject takes some steps to rotate back to the sagittal plane (idxF3:idxF4)
    % 5: Subject deposits the box in the sagittal plane (idxF4:idxF5)
    

%% Duration of each phase

% phase 1 
phase1 = LTOE_z(idxO1:idxF1);
t_phase1= double(length(phase1))*Ts;

% phase 2
phase2 = LTOE_z(idxF1:idxF2);
t_phase2= double(length(phase2)-1)*Ts;

% phase 3
phase3 = LTOE_z(idxF2:idxF3);
t_phase3= double(length(phase3)-1)*Ts;

% phase 4
phase4 = LTOE_z(idxF3:idxF4);
t_phase4= double(length(phase4)-1)*Ts;

% phase 5
phase5 = LTOE_z(idxF4:idxF5);
t_phase5= double(length(phase5)-1)*Ts;


%% Verification
Sum = t_phase1 + t_phase2 + t_phase3 + t_phase4 + t_phase5;

if round(Sum - t_total) == 0
    disp('Está bien hecho!');
else 
    disp('Algo ha fallado');
end

Time = table([t_phase1; t_phase2; t_phase3; t_phase4; t_phase5; Sum; t_total],'VariableNames',{'Seconds'},'RowNames',{'Phase 1','Phase 2','Phase 3', 'Phase 4', 'Phase 5', 'Sum', 'Complete Trial'})
