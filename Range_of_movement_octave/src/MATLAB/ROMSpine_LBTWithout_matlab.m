%% Range of movement of Spine in Lateral Box Transfer

% Code AUTHOR: Yaiza Benito Molpeceres. DATE: January-May 2020.

% Without exoskeleton trials segmentation: din 42, 43, 47, 48, 49, 08, 21.
% Three phases shown in figure 1:
    % 1: Subject picks up the box in the sagittal plane and takes a step to rotate to the frontal plane (idxO1:idxF2)
    % 2: Subject deposits the box in the frontal plane and takes a step to rotate back to the sagittal plane(idxF2:idxF4)
    % 3: Subject deposits the box in the sagittal plane (idxF4:idxF5)

clear all % Clear variables
close all % Close figures
clc

load('..\tests\data\input\dinamica42_B.mat')
             
Ts = 1/Fs;
t_total = (double(frames)*Ts);
t_1000 = 0:(t_total/1000):(t_total-(t_total/1000));

% Normal Ranges of Spine Motion
for i =1:(length(t_1000))
    flex_spine(1,i) = 75; % flexion
    ext_spine(1,i) = -30;   % extension
    abd_spine(1,i) = -35;   % abduction
end

% Feet markers for signal segmentation
%LTOE_z = LTOE(:,3)'; 
%RTOE_z = RTOE(:,3)';   
LHEE_z = LHEE(:,3)';  
RHEE_z = RHEE(:,3)';    
  

%% SEGMENTATION FIRST TURN: 

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
xlabel('Time in frames');
ylabel('Milimetres');
title('Segmentation Without Exoskeleton Trials')

%% Import model data for Spine Relative Angles

% Spine  
        RSpineAngles_x = ModelData.Raw.(ModelOutput{11})(1,:); 
        RSpineAngles_y = ModelData.Raw.(ModelOutput{11})(2,:);
        RSpineAngles_z = ModelData.Raw.(ModelOutput{11})(3,:);
        

        LSpineAngles_x = ModelData.Raw.(ModelOutput{12})(1,:);
        LSpineAngles_y = ModelData.Raw.(ModelOutput{12})(2,:);
        LSpineAngles_z = ModelData.Raw.(ModelOutput{12})(3,:);
        

%% Spine signals visualization before going through segmentation

figure(2)
frames = length(LHEE_z);
% (1:(frames/1000):frames)  si frames>1000
% (1:(frames/1000):frames+(frames/1000)) si frames<1000
% +(frames/1000)
plot(t_1000,RSpineAngles_x(1:(frames/1000):frames)); 
hold on 
plot(t_1000, flex_spine,'--',t_1000, ext_spine,'--'); 
legend({'x', 'Normal range of flexion-extension'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
%ylim([-100 160]);
title('Right Spine angles: Flexo-Extension');

figure(3)
plot(t_1000, RSpineAngles_y(1:(frames/1000):frames));
hold on 
plot(t_1000, abd_spine,'--'); 
legend({'y', 'Normal range of abduction'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
ylim([-40 10]);
title('Right Spine angles: Abduction-Adduction');

figure(4)
plot(t_1000, RSpineAngles_z(1:(frames/1000):frames));
legend({'z'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
%ylim([-70 50]);
title('Right Spine angles: Internal-External Rotation');

figure(5)
% (1:(frames/1000):frames)  si frames>1000
% (1:(frames/1000):frames+(frames/1000)) si frames<1000
% +(frames/1000)
plot(t_1000,LSpineAngles_x(1:(frames/1000):frames)); 
hold on 
plot(t_1000, flex_spine,'--',t_1000, ext_spine,'--'); 
legend({'x', 'Normal range of flexion-extension'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
%ylim([-100 160]);
title('Left Spine angles: Flexo-Extension');

figure(6)
plot(t_1000, LSpineAngles_y(1:(frames/1000):frames));
hold on 
plot(t_1000, abd_spine,'--'); 
legend({'y', 'Normal range of abduction'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
ylim([-40 5]);
title('Left Spine angles: Abduction-Adduction');

figure(7)
plot(t_1000, LSpineAngles_z(1:(frames/1000):frames));
legend({'z'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
%ylim([-70 80]);
title('Left Spine angles: Internal-External Rotation');

%figure(3)
%plot(t ,RSpineAngles_x, t, RSpineAngles_y, t, RSpineAngles_z);
%legend('x','y','z');
%xlabel('time in frames'); % frames totales del trial
%title('Right Spine angles');

%figure(4)
%plot(t,LSpineAngles_x, t, LSpineAngles_y, t, LSpineAngles_z);   
%legend('x','y','z');
%xlabel('time in frames'); % frames totales del trial
%title('Left Spine angles');

%% Signal Segmentation

% Right Spine 

phase1_RSpine_x = RSpineAngles_x(idxO1:idxF2);  
phase1_RSpine_y = RSpineAngles_y(idxO1:idxF2);
phase1_RSpine_z = RSpineAngles_z(idxO1:idxF2);

phase2_RSpine_x = RSpineAngles_x(idxF2:idxF4); 
phase2_RSpine_y = RSpineAngles_y(idxF2:idxF4);
phase2_RSpine_z = RSpineAngles_z(idxF2:idxF4);

phase3_RSpine_x = RSpineAngles_x(idxF4:idxF5); 
phase3_RSpine_y = RSpineAngles_y(idxF4:idxF5);
phase3_RSpine_z = RSpineAngles_z(idxF4:idxF5);

% Left Spine

phase1_LSpine_x = LSpineAngles_x(idxO1:idxF2);
phase1_LSpine_y = LSpineAngles_y(idxO1:idxF2);
phase1_LSpine_z = LSpineAngles_z(idxO1:idxF2);

phase2_LSpine_x = LSpineAngles_x(idxF2:idxF4);
phase2_LSpine_y = LSpineAngles_y(idxF2:idxF4);
phase2_LSpine_z = LSpineAngles_z(idxF2:idxF4);

phase3_LSpine_x = LSpineAngles_x(idxF4:idxF5);
phase3_LSpine_y = LSpineAngles_y(idxF4:idxF5);
phase3_LSpine_z = LSpineAngles_z(idxF4:idxF5);

%% time vectors definition

t2 = 0:1:(length(phase3_RSpine_x)-1);
t3 = 0:1:(length(phase1_RSpine_x)-1);
t4 = 0:1:(length(phase2_RSpine_x)-1);

% Right Spine figures: MOVEMENT SEGMENTED INTO THE 4 PHASES

figure(8)
subplot(3,1,1)
plot (t3, phase1_RSpine_x, t3, phase1_RSpine_y, '--',t3,phase1_RSpine_z, '.');
legend('x','y','z')
xlabel('Time in frames');  
ylabel('Degrees');
title('RIGHT ANGLES: Phase 1'); % Box Carrying from Sagittal to Frontal
subplot(3,1,2)
plot(t4, phase2_RSpine_x,t4,phase2_RSpine_y, '--',t4,phase2_RSpine_z, '.');
legend('x','y','z')
xlabel('Time in frames');  
ylabel('Degrees');
title('Phase 2'); % Box Carrying from Frontal to Sagittal
subplot(3,1,3)
plot(t2, phase3_RSpine_x, t2, phase3_RSpine_y,'--', t2, phase3_RSpine_z, '.');
legend('x','y','z')
xlabel('Time in frames');  
ylabel('Degrees');
title('Phase 3'); % Placing box sagittal plane again


% Left Spine figures: MOVEMENT SEGMENTED INTO THE 4 PHASES

figure(9)
subplot(3,1,1)
plot (t3, phase1_LSpine_x, t3, phase1_LSpine_y, '--',t3,phase1_LSpine_z, '.');
legend('x','y','z')
xlabel('Time in frames');  
ylabel('Degrees');
title('LEFT ANGLES: Phase 1'); % Box Carrying from Sagittal to Frontal
subplot(3,1,2)
plot(t4, phase2_LSpine_x,t4, phase2_LSpine_y, '--',t4, phase2_LSpine_z, '.');
legend('x','y','z')
xlabel('Time in frames');  
ylabel('Degrees');
title('Phase 2'); % Box Carrying from Frontal to Sagittal
subplot(3,1,3)
plot(t2, phase3_LSpine_x, t2, phase3_LSpine_y,'--', t2, phase3_LSpine_z, '.');
legend('x','y','z')
xlabel('Time in frames');  
ylabel('Degrees');
title('Phase 3'); % Placing box sagittal plane again


%% ROM calculation

% Right Spine

RSpine_ph1(:,1) = phase1_RSpine_x;
RSpine_ph1(:,2) = phase1_RSpine_y;
RSpine_ph1(:,3) = phase1_RSpine_z;

RSpine_ph2(:,1) = phase2_RSpine_x;
RSpine_ph2(:,2) = phase2_RSpine_y;
RSpine_ph2(:,3) = phase2_RSpine_z;

RSpine_ph3(:,1) = phase3_RSpine_x;
RSpine_ph3(:,2) = phase3_RSpine_y;
RSpine_ph3(:,3) = phase3_RSpine_z;


RSpine_ph1_max = max(RSpine_ph1);
RSpine_ph1_min = min(RSpine_ph1);

RSpine_ph2_max = max(RSpine_ph2);
RSpine_ph2_min = min(RSpine_ph2);

RSpine_ph3_max = max(RSpine_ph3);
RSpine_ph3_min = min(RSpine_ph3);


for i = 1:length(RSpine_ph1_max)
    RSpine_ph1_ROM(1,i) = RSpine_ph1_max(i)- RSpine_ph1_min(i);
end

for i = 1:length(RSpine_ph1_max)
    RSpine_ph2_ROM(1,i) = RSpine_ph2_max(i)- RSpine_ph2_min(i);
end

for i = 1:length(RSpine_ph1_max)
    RSpine_ph3_ROM(1,i) = RSpine_ph3_max(i)- RSpine_ph3_min(i);
end

% Left Spine

LSpine_ph1(:,1) = phase1_LSpine_x;
LSpine_ph1(:,2) = phase1_LSpine_y;
LSpine_ph1(:,3) = phase1_LSpine_z;

LSpine_ph2(:,1) = phase2_LSpine_x;
LSpine_ph2(:,2) = phase2_LSpine_y;
LSpine_ph2(:,3) = phase2_LSpine_z;

LSpine_ph3(:,1) = phase3_LSpine_x;
LSpine_ph3(:,2) = phase3_LSpine_y;
LSpine_ph3(:,3) = phase3_LSpine_z;


LSpine_ph1_max = max(LSpine_ph1);
LSpine_ph1_min = min(LSpine_ph1);

LSpine_ph2_max = max(LSpine_ph2);
LSpine_ph2_min = min(LSpine_ph2);

LSpine_ph3_max = max(LSpine_ph3);
LSpine_ph3_min = min(LSpine_ph3);


for i = 1:length(RSpine_ph1_max)
    LSpine_ph1_ROM(1,i) = LSpine_ph1_max(i)- LSpine_ph1_min(i);
end

for i = 1:length(RSpine_ph1_max)
    LSpine_ph2_ROM(1,i) = LSpine_ph2_max(i)- LSpine_ph2_min(i);
end

for i = 1:length(RSpine_ph1_max)
    LSpine_ph3_ROM(1,i) = LSpine_ph3_max(i)- LSpine_ph3_min(i);
end


% Total ROM Spine Angles

RSpine_total(:,1) = RSpineAngles_x;
RSpine_total(:,2) = RSpineAngles_y;
RSpine_total(:,3) = RSpineAngles_z;

RSpine_total_max = max(RSpine_total);
RSpine_total_min = min(RSpine_total);

for i = 1:length(RSpine_total_max)
    RSpine_total_ROM(1,i) = RSpine_total_max(i)- RSpine_total_min(i);
end

LSpine_total(:,1) = LSpineAngles_x;
LSpine_total(:,2) = LSpineAngles_y;
LSpine_total(:,3) = LSpineAngles_z;

LSpine_total_max = max(LSpine_total);
LSpine_total_min = min(LSpine_total);

for i = 1:length(LSpine_total_max)
    LSpine_total_ROM(1,i) = LSpine_total_max(i)- LSpine_total_min(i);
end

ROM_Right_Spine = table([RSpine_ph1_ROM(1);RSpine_ph2_ROM(1);RSpine_ph3_ROM(1)],[RSpine_ph1_ROM(2);RSpine_ph2_ROM(2);RSpine_ph3_ROM(2)],[RSpine_ph1_ROM(3);RSpine_ph2_ROM(3);RSpine_ph3_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Phase 1','Phase 2','Phase 3'})

ROM_Left_Spine = table([LSpine_ph1_ROM(1);LSpine_ph2_ROM(1);LSpine_ph3_ROM(1)],[LSpine_ph1_ROM(2);LSpine_ph2_ROM(2);LSpine_ph3_ROM(2)],[LSpine_ph1_ROM(3);LSpine_ph2_ROM(3);LSpine_ph3_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Phase 1','Phase 2','Phase 3'})

TotalROM_Spine = table([RSpine_total_ROM(1); LSpine_total_ROM(1)],[RSpine_total_ROM(2); LSpine_total_ROM(2)],[RSpine_total_ROM(3); LSpine_total_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Right','Left'})

