%% Total time to perform the task in Lateral Box Transfer (without exoskeleton segmentation)

% Code AUTHOR: Yaiza Benito Molpeceres. DATE: January-May 2020.

clear all % Clear variables
close all % Close figures
clc

load('../tests/data/input/dinamica42_B.mat')
             
Ts = 1/Fs;

%% Total time to perform the task 

t_total = (double(frames)*Ts);

%% Events of interest in feet markers trajectories

% z LHEE RHEE marker trajectories
LHEE_z = LHEE(:,3)';  
RHEE_z = RHEE(:,3)'; 

% First turn

[pks, locs] = findpeaks(LHEE_z,'minPeakProminence',10);
TF1 = islocalmin(LHEE_z, 'FlatSelection', 'first');
idx = find(TF1);
flat2 = idx > locs(1);
idx_flat2 = find(flat2);

[pks2, locs2] = findpeaks(RHEE_z,'minPeakProminence',10);
TF2 = islocalmin(RHEE_z, 'FlatSelection', 'first');
idx2 = find(TF2);
flat3 = idx2 < locs2(1);  % 42 es <
idx_flat3 = find(flat3);

% Second turn

flat5 = idx < locs(2);      
idx_flat5 = find(flat5);

flat5 = idx > locs(2);
idx_flat5 = find(flat5);

flat8 = idx2 > locs2(2);
idx_flat8 = find(flat8);

% Indexes definition for signal segmentation
idxO1 = 1;
idxF1 = idx2(idx_flat3(length(idx_flat3))); 
idxF2 = idx(idx_flat2(1)); 
idxF3 = idx(idx_flat5(length(idx_flat5))); 
idxF4 = idx2(idx_flat8(1)); 
idxF5 = length(LHEE_z);
   
% Visualization of the segmentation

figure(1)
t= 0:(length(LHEE)-1);
plot(t, zscore(LHEE_z));
hold on
plot(t,zscore(RHEE_z));
plot(t,zscore(LHEE_z),'r*', 'MarkerIndices', idxO1);
plot(t,zscore(LHEE_z),'r*', 'MarkerIndices', idxF2); % 42 idxF2
plot(t,zscore(LHEE_z),'r*', 'MarkerIndices', idxF5);
plot(t,zscore(RHEE_z),'r*', 'MarkerIndices', idxO1);
plot(t,zscore(RHEE_z),'r*', 'MarkerIndices', idxF4); % 42 idxF4
plot(t,zscore(RHEE_z),'r*', 'MarkerIndices', idxF5);
hold off
legend('Left Heel Z', 'Right Heel Z');
title('Segmentation Without Exoskeleton Trials')

% Three phases:
    % 1: Subject picks up the box in the sagittal plane and takes a step to rotate to the frontal plane (idxO1:idxF2)
    % 2: Subject deposits the box in the frontal plane and takes a step to rotate back to the sagittal plane(idxF2:idxF4)
    % 3: Subject deposits the box in the sagittal plane (idxF4:idxF5)
    

%% Duration of each phase

% phase 1 
phase1 = LHEE_z(idxO1:idxF2);
t_phase1= double(length(phase1))*Ts;

% phase 2
phase2 = LHEE_z(idxF2:idxF4);
t_phase2= double(length(phase2)-1)*Ts;

% phase 3
phase3 = LHEE_z(idxF4:idxF5);
t_phase3= double(length(phase3)-1)*Ts;


%% Verification
Sum = t_phase1 + t_phase2 + t_phase3;

if round(Sum - t_total) == 0
    disp('Está bien hecho!');
else 
    disp('Algo ha fallado');
end

Time = table([t_phase1; t_phase2; t_phase3; Sum; t_total],'VariableNames',{'Seconds'},'RowNames',{'Phase 1','Phase 2','Phase 3', 'Sum', 'Complete Trial'})




