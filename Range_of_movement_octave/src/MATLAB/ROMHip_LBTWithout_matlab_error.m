%% Range of movement of Hip in Lateral Box Transfer

% Code AUTHOR: Yaiza Benito Molpeceres. DATE: January-May 2020.

% Without exoskeleton trials segmentation: din 42, 43, 47, 48, 49, 08, 21.
% Three phases shown in figure 1:
    % 1: Subject picks up the box in the sagittal plane and takes a step to rotate to the frontal plane (idxO1:idxF2),(idxO1:idxF1)
    % 2: Subject deposits the box in the frontal plane and takes a step to rotate back to the sagittal plane (idxF2:idxF4),(idxF1:idxF3)
    % 3: Subject deposits the box in the sagittal plane (idxF4:idxF5),(idxF3:idxF5)

clear all % Clear variables
close all % Close figures
clc

load('..\tests\data\input\dinamica08_B.mat')
             
Ts = 1/Fs;
t_total = (double(frames)*Ts);
t_1000 = 0:(t_total/1000):(t_total-(t_total/1000));

% Normal Ranges of Hip Motion
for i =1:(length(t_1000))
    flex_hip(1,i) = 130; % flexion
    ext_hip(1,i) = -30;   % extension
    abd_hip(1,i) = -50;   % abduction
    add_hip(1,i) = 30;   % adduction
    intRot_hip(1,i) = 40;   % internal rotation
    extRot_hip(1,i) = -45;   % external rotation
end

% Feet markers for signal segmentation
%LTOE_z = LTOE(:,3)'; 
%RTOE_z = RTOE(:,3)';   
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
flat3 = idx2 > locs2(1);  
idx_flat3 = find(flat3);

% Second turn

%flat5 = idx < locs(2);      
%idx_flat5 = find(flat5);

flat5 = idx > locs(3);
idx_flat5 = find(flat5);

flat8 = idx2 > locs2(2);
idx_flat8 = find(flat8);

% Indexes definition for signal segmentation
idxO1 = 1;
idxF1 = idx2(idx_flat3(1)); % pie derecho 
idxF2 = idx(idx_flat2(1)); 
idxF3 = idx(idx_flat5(1)); % pie izquierdo
idxF4 = idx2(idx_flat8(1)); 
idxF5 = length(LHEE_z);
   
% Visualization of the segmentation

figure(1)
t= 0:(length(LHEE)-1);
plot(t, zscore(LHEE_z));
hold on
plot(t,zscore(RHEE_z));
plot(t,zscore(LHEE_z),'r*', 'MarkerIndices', idxO1);
plot(t,zscore(LHEE_z),'r*', 'MarkerIndices', idxF3); 
plot(t,zscore(LHEE_z),'r*', 'MarkerIndices', idxF5);
plot(t,zscore(RHEE_z),'r*', 'MarkerIndices', idxO1);
plot(t,zscore(RHEE_z),'r*', 'MarkerIndices', idxF1); 
plot(t,zscore(RHEE_z),'r*', 'MarkerIndices', idxF5);
hold off
legend('Left Heel Z', 'Right Heel Z');
title('Segmentation Without Exoskeleton Trials')

%% Import model data for Hip Relative Angles

% Hip  
        RHipAngles_x = ModelData.Raw.(ModelOutput{1})(1,:); 
        RHipAngles_y = ModelData.Raw.(ModelOutput{1})(2,:);
        RHipAngles_z = ModelData.Raw.(ModelOutput{1})(3,:);
        

        LHipAngles_x = ModelData.Raw.(ModelOutput{2})(1,:);
        LHipAngles_y = ModelData.Raw.(ModelOutput{2})(2,:);
        LHipAngles_z = ModelData.Raw.(ModelOutput{2})(3,:);
        

%% Hip signals visualization before going through segmentation

figure(2)
frames = length(LHEE_z);
% (1:(frames/1000):frames)  si frames>1000
% (1:(frames/1000):frames+(frames/1000)) si frames<1000
% +(frames/1000)
plot(t_1000,RHipAngles_x(1:(frames/1000):frames)); 
hold on 
plot(t_1000, flex_hip,'--',t_1000, ext_hip,'--'); 
legend({'x', 'Normal range of flexion-extension'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
%ylim([-100 160]);
title('Right Hip angles: Flexo-Extension');

figure(3)
plot(t_1000, RHipAngles_y(1:(frames/1000):frames));
hold on 
plot(t_1000, abd_hip,'--',t_1000, add_hip,'--'); 
legend({'y', 'Normal range of abduction-adduction'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
ylim([-60 40]);
title('Right Hip angles: Abduction-Adduction');

figure(4)
plot(t_1000, RHipAngles_z(1:(frames/1000):frames));
hold on 
plot(t_1000, intRot_hip,'--',t_1000, extRot_hip,'--'); 
legend({'z', 'Normal range of internal-external rotation'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
ylim([-50 50]);
title('Right Hip angles: Internal-External Rotation');

figure(5)
% (1:(frames/1000):frames)  si frames>1000
% (1:(frames/1000):frames+(frames/1000)) si frames<1000
% +(frames/1000)
plot(t_1000,LHipAngles_x(1:(frames/1000):frames));
hold on 
plot(t_1000, flex_hip,'--',t_1000, ext_hip,'--'); 
legend({'x', 'Normal range of flexion-extension'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
%ylim([-100 160]);
title('Left Hip angles: Flexo-Extension');

figure(6)
plot(t_1000, LHipAngles_y(1:(frames/1000):frames));
hold on 
plot(t_1000, abd_hip,'--',t_1000, add_hip,'--'); 
legend({'y', 'Normal range of abduction-adduction'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
ylim([-60 40]);
title('Left Hip angles: Abduction-Adduction');

figure(7)
plot(t_1000, LHipAngles_z(1:(frames/1000):frames));
hold on 
plot(t_1000, intRot_hip,'--',t_1000, extRot_hip,'--'); 
legend({'z', 'Normal range of internal-external rotation'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
%ylim([-70 80]);
title('Left Hip angles: Internal-External Rotation');

%figure(3)
%plot(t ,RHipAngles_x, t, RHipAngles_y, t, RHipAngles_z);
%legend('x','y','z');
%xlabel('time in frames'); % frames totales del trial
%title('Right Hip angles');

%figure(4)
%plot(t,LHipAngles_x, t, LHipAngles_y, t, LHipAngles_z);   
%legend('x','y','z');
%xlabel('time in frames'); % frames totales del trial
%title('Left Hip angles');

%% Signal Segmentation

% Right Hip 

phase1_RHip_x = RHipAngles_x(idxO1:idxF1);  
phase1_RHip_y = RHipAngles_y(idxO1:idxF1);
phase1_RHip_z = RHipAngles_z(idxO1:idxF1);

phase2_RHip_x = RHipAngles_x(idxF1:idxF3); 
phase2_RHip_y = RHipAngles_y(idxF1:idxF3);
phase2_RHip_z = RHipAngles_z(idxF1:idxF3);

phase3_RHip_x = RHipAngles_x(idxF3:idxF5); 
phase3_RHip_y = RHipAngles_y(idxF3:idxF5);
phase3_RHip_z = RHipAngles_z(idxF3:idxF5);

% Left Hip

phase1_LHip_x = LHipAngles_x(idxO1:idxF1);
phase1_LHip_y = LHipAngles_y(idxO1:idxF1);
phase1_LHip_z = LHipAngles_z(idxO1:idxF1);

phase2_LHip_x = LHipAngles_x(idxF1:idxF3);
phase2_LHip_y = LHipAngles_y(idxF1:idxF3);
phase2_LHip_z = LHipAngles_z(idxF1:idxF3);

phase3_LHip_x = LHipAngles_x(idxF3:idxF5);
phase3_LHip_y = LHipAngles_y(idxF3:idxF5);
phase3_LHip_z = LHipAngles_z(idxF3:idxF5);

%% time vectors definition

t2 = 0:1:(length(phase3_RHip_x)-1);
t3 = 0:1:(length(phase1_RHip_x)-1);
t4 = 0:1:(length(phase2_RHip_x)-1);

% Right Hip figures: MOVEMENT SEGMENTED INTO THE 4 PHASES

figure(8)
subplot(3,1,1)
plot (t3, phase1_RHip_x, t3, phase1_RHip_y, '--',t3,phase1_RHip_z, '.');
legend('x','y','z')
xlabel('Time in frames'); 
ylabel('Degrees');
title('RIGHT HIP: Phase 1');   % Box Carrying from Sagittal to Frontal
subplot(3,1,2)
plot(t4, phase2_RHip_x,t4,phase2_RHip_y, '--',t4,phase2_RHip_z, '.');
legend('x','y','z')
xlabel('Time in frames');  
ylabel('Degrees');
title('Phase 2');  % Box Carrying from Frontal to Sagittal
subplot(3,1,3)
plot(t2, phase3_RHip_x, t2, phase3_RHip_y,'--', t2, phase3_RHip_z, '.');
legend('x','y','z')
xlabel('Time in frames'); 
ylabel('Degrees');
title('Phase 3');  % Placing box sagittal plane again


% Left Hip figures: MOVEMENT SEGMENTED INTO THE 4 PHASES

figure(9)
subplot(3,1,1)
plot (t3, phase1_LHip_x, t3, phase1_LHip_y, '--',t3,phase1_LHip_z, '.');
legend('x','y','z')
xlabel('Time in frames'); 
ylabel('Degrees');
title('LEFT HIP: Phase 1');    % Box Carrying from Sagittal to Frontal
subplot(3,1,2)
plot(t4, phase2_LHip_x,t4, phase2_LHip_y, '--',t4, phase2_LHip_z, '.');
legend('x','y','z')
xlabel('Time in frames');  
ylabel('Degrees');
title('Phase 2');  % Box Carrying from Frontal to Sagittal
subplot(3,1,3)
plot(t2, phase3_LHip_x, t2, phase3_LHip_y,'--', t2, phase3_LHip_z, '.');
legend('x','y','z')
xlabel('Time in frames'); 
ylabel('Degrees');
title('Phase 3');  % Placing box sagittal plane again


%% ROM calculation

% Right Hip

RHip_ph1(:,1) = phase1_RHip_x;
RHip_ph1(:,2) = phase1_RHip_y;
RHip_ph1(:,3) = phase1_RHip_z;

RHip_ph2(:,1) = phase2_RHip_x;
RHip_ph2(:,2) = phase2_RHip_y;
RHip_ph2(:,3) = phase2_RHip_z;

RHip_ph3(:,1) = phase3_RHip_x;
RHip_ph3(:,2) = phase3_RHip_y;
RHip_ph3(:,3) = phase3_RHip_z;


RHip_ph1_max = max(RHip_ph1);
RHip_ph1_min = min(RHip_ph1);

RHip_ph2_max = max(RHip_ph2);
RHip_ph2_min = min(RHip_ph2);

RHip_ph3_max = max(RHip_ph3);
RHip_ph3_min = min(RHip_ph3);


for i = 1:length(RHip_ph1_max)
    RHip_ph1_ROM(1,i) = RHip_ph1_max(i)- RHip_ph1_min(i);
end

for i = 1:length(RHip_ph1_max)
    RHip_ph2_ROM(1,i) = RHip_ph2_max(i)- RHip_ph2_min(i);
end

for i = 1:length(RHip_ph1_max)
    RHip_ph3_ROM(1,i) = RHip_ph3_max(i)- RHip_ph3_min(i);
end

% Left Hip

LHip_ph1(:,1) = phase1_LHip_x;
LHip_ph1(:,2) = phase1_LHip_y;
LHip_ph1(:,3) = phase1_LHip_z;

LHip_ph2(:,1) = phase2_LHip_x;
LHip_ph2(:,2) = phase2_LHip_y;
LHip_ph2(:,3) = phase2_LHip_z;

LHip_ph3(:,1) = phase3_LHip_x;
LHip_ph3(:,2) = phase3_LHip_y;
LHip_ph3(:,3) = phase3_LHip_z;


LHip_ph1_max = max(LHip_ph1);
LHip_ph1_min = min(LHip_ph1);

LHip_ph2_max = max(LHip_ph2);
LHip_ph2_min = min(LHip_ph2);

LHip_ph3_max = max(LHip_ph3);
LHip_ph3_min = min(LHip_ph3);


for i = 1:length(RHip_ph1_max)
    LHip_ph1_ROM(1,i) = LHip_ph1_max(i)- LHip_ph1_min(i);
end

for i = 1:length(RHip_ph1_max)
    LHip_ph2_ROM(1,i) = LHip_ph2_max(i)- LHip_ph2_min(i);
end

for i = 1:length(RHip_ph1_max)
    LHip_ph3_ROM(1,i) = LHip_ph3_max(i)- LHip_ph3_min(i);
end


% Total ROM Hip Angles

RHip_total(:,1) = RHipAngles_x;
RHip_total(:,2) = RHipAngles_y;
RHip_total(:,3) = RHipAngles_z;

RHip_total_max = max(RHip_total);
RHip_total_min = min(RHip_total);

for i = 1:length(RHip_total_max)
    RHip_total_ROM(1,i) = RHip_total_max(i)- RHip_total_min(i);
end

LHip_total(:,1) = LHipAngles_x;
LHip_total(:,2) = LHipAngles_y;
LHip_total(:,3) = LHipAngles_z;

LHip_total_max = max(LHip_total);
LHip_total_min = min(LHip_total);

for i = 1:length(LHip_total_max)
    LHip_total_ROM(1,i) = LHip_total_max(i)- LHip_total_min(i);
end

ROM_Right_Hip = table([RHip_ph1_ROM(1);RHip_ph2_ROM(1);RHip_ph3_ROM(1)],[RHip_ph1_ROM(2);RHip_ph2_ROM(2);RHip_ph3_ROM(2)],[RHip_ph1_ROM(3);RHip_ph2_ROM(3);RHip_ph3_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Phase 1','Phase 2','Phase 3'})

ROM_Left_Hip = table([LHip_ph1_ROM(1);LHip_ph2_ROM(1);LHip_ph3_ROM(1)],[LHip_ph1_ROM(2);LHip_ph2_ROM(2);LHip_ph3_ROM(2)],[LHip_ph1_ROM(3);LHip_ph2_ROM(3);LHip_ph3_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Phase 1','Phase 2','Phase 3'})

TotalROM_Hip = table([RHip_total_ROM(1); LHip_total_ROM(1)],[RHip_total_ROM(2); LHip_total_ROM(2)],[RHip_total_ROM(3); LHip_total_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Right','Left'})

