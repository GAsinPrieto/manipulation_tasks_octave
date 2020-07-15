%% Range of movement of Elbow in Lateral Box Transfer

% Code AUTHOR: Yaiza Benito Molpeceres. DATE: January-May 2020.

% Exoskeleton trials segmentation: 54, 55, 56, 62.
% Five phases shown in figure 1:
    % 1: Subject picks up the box in the sagittal plane (idxO1:idxF1)
    % 2: Subject takes some steps to rotate to the frontal plane (idxF1:idxF2)
    % 3: Subject deposits the box in the frontal plane (idxF2:idxF3)
    % 4: Subject takes some steps to rotate back to the sagittal plane (idxF3:idxF4)
    % 5: Subject deposits the box in the sagittal plane (idxF4:idxF5)


clear all % Clear variables
close all % Close figures
clc

load('..\tests\data\input\dinamica56_B.mat')
           
Ts = 1/Fs;
t_total = (double(frames)*Ts);
t_1000 = 0:(t_total/1000):(t_total-(t_total/1000));

% Normal Ranges of Elbow Motion
for i =1:(length(t_1000))
    flex_elbow(1,i) = 130; % flexion o 180 V
    ext_elbow(1,i) = -45;   % extension o 60 V
    supination_elbow(1,i) = 20;  % supination
    pronation_elbow(1,i) = -20;  % pronation
end

% Feet markers for signal segmentation
LTOE_z = LTOE(:,3)'; 
RTOE_z = RTOE(:,3)';   
%LHEE_z = LHEE(:,3)';  
%RHEE_z = RHEE(:,3)'; 

t= 0:(length(LTOE)-1);

%% SEGMENTATION LEFT TOE: 

[pks, locs] = findpeaks(LTOE_z,'minPeakProminence',10,'MinPeakHeight',80);                     
TF1 = islocalmin(LTOE_z, 'FlatSelection', 'first');
idx = find(TF1);
flat = idx < locs(1);
idx_flat = find(flat);

flat2 = idx < locs(3)& idx > locs(2);
idx_flat2 = find(flat2);

flat3 = idx > locs(4);
idx_flat3 = find(flat3);


figure(1)
plot(t, LTOE_z, t(locs), LTOE_z(locs),'o')
hold on
plot(t, LTOE_z,'r*','MarkerIndices',idx(idx_flat(length(idx_flat))));
plot(t, LTOE_z,'r*','MarkerIndices',idx(idx_flat2(3)));                    
plot(t, LTOE_z,'r*','MarkerIndices',idx(idx_flat2(length(idx_flat2))));
plot(t, LTOE_z,'r*','MarkerIndices',idx(idx_flat3(4)));                     
hold off
title('Left Toe');

%% SEGMENTATION RIGHT TOE: 

[pks2, locs2] = findpeaks(RTOE_z,'minPeakProminence',10);

TF2 = islocalmin(RTOE_z, 'FlatSelection', 'first');
idx2 = find(TF2);
flat4 = idx2 < locs2(1);
idx_flat4 = find(flat4);

flat5 = idx2 > locs2(length(locs2));
idx_flat5 = find(flat5);
flat6 = idx2 > locs2(3);
idx_flat6 = find(flat6);

figure(2)
plot(t, RTOE_z, t(locs2), RTOE_z(locs2),'o')
hold on
plot(t, RTOE_z,'r*','MarkerIndices',idx2(idx_flat4(length(idx_flat4))));
plot(t, RTOE_z,'r*','MarkerIndices',idx2(idx_flat6(1)));
plot(t, RTOE_z,'r*','MarkerIndices',idx2(idx_flat5(2))); 
title('Right Toe');

% Indexes definiton for signal segmentation
idxO1 = 1;
idxF1 = idx2(idx_flat4(length(idx_flat4))); %idxO2
idxF2 = idx2(idx_flat6(1)); %idxO3
idxF3 = idx(idx_flat2(length(idx_flat2))); %idxO4
idxF4 = idx(idx_flat3(3)); %idxO5
idxF5 = length(LTOE_z);

figure(3)
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
xlabel('Time in frames');
ylabel('Milimetres');
title('Lateral Box Transfer Segmentation for Exoskeleton Trials')


%% Import model data for Elbow Relative Angles

% Elbow  
        RElbowAngles_x = ModelData.Raw.(ModelOutput{7})(1,:); 
        RElbowAngles_y = ModelData.Raw.(ModelOutput{7})(2,:);
        RElbowAngles_z = ModelData.Raw.(ModelOutput{7})(3,:);
        

        LElbowAngles_x = ModelData.Raw.(ModelOutput{8})(1,:);
        LElbowAngles_y = ModelData.Raw.(ModelOutput{8})(2,:);
        LElbowAngles_z = ModelData.Raw.(ModelOutput{8})(3,:);
        

%% Elbow signals visualization before going through segmentation

figure(4)
frames = length(LTOE_z);
% (1:(frames/1000):frames)  si frames>1000
% (1:(frames/1000):frames+(frames/1000)) si frames<1000
% +(frames/1000)
plot(t_1000,RElbowAngles_x(1:(frames/1000):frames)); 
hold on 
plot(t_1000, flex_elbow,'--',t_1000, ext_elbow,'--'); 
legend({'x', 'Normal range of flexion-extension'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
%ylim([-100 160]);
title('Right Elbow angles: Flexo-Extension');

figure(5)
plot(t_1000, RElbowAngles_y(1:(frames/1000):frames));
legend({'y'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
%ylim([-70 50]);
title('Right Elbow angles: Supination-Pronation');

figure(6)
plot(t_1000, RElbowAngles_z(1:(frames/1000):frames));
legend({'z'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
%ylim([-70 50]);
title('Right Elbow angles: Internal-External Rotation');

figure(7)
% (1:(frames/1000):frames)  si frames>1000
% (1:(frames/1000):frames+(frames/1000)) si frames<1000
% +(frames/1000)
plot(t_1000,LElbowAngles_x(1:(frames/1000):frames)); 
hold on 
plot(t_1000, flex_elbow,'--',t_1000, ext_elbow,'--'); 
legend({'x', 'Normal range of flexion-extension'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
%ylim([-100 160]);
title('Left Elbow angles: Flexo-Extension');

figure(8)
plot(t_1000, LElbowAngles_y(1:(frames/1000):frames));
legend({'y'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
%ylim([-70 50]);
title('Left Elbow angles: Supination-Pronation');

figure(9)
plot(t_1000, LElbowAngles_z(1:(frames/1000):frames));
legend({'z'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
%ylim([-70 80]);
title('Left Elbow angles: Internal-External Rotation');


%figure(3)
%plot(t ,RElbowAngles_x, t, RElbowAngles_y, t, RElbowAngles_z);
%legend('x','y','z');
%xlabel('time in frames'); % frames totales del trial
%title('Right Elbow angles');

%figure(4)
%plot(t,LElbowAngles_x, t, LElbowAngles_y, t, LElbowAngles_z);   
%legend('x','y','z');
%xlabel('time in frames'); % frames totales del trial
%title('Left Elbow angles');

%% Signal Segmentation

% Right Elbow 

phase1_RElbow_x = RElbowAngles_x(idxO1:idxF1);  
phase1_RElbow_y = RElbowAngles_y(idxO1:idxF1);
phase1_RElbow_z = RElbowAngles_z(idxO1:idxF1);

phase2_RElbow_x = RElbowAngles_x(idxF1:idxF2); 
phase2_RElbow_y = RElbowAngles_y(idxF1:idxF2);
phase2_RElbow_z = RElbowAngles_z(idxF1:idxF2);

phase3_RElbow_x = RElbowAngles_x(idxF2:idxF3);
phase3_RElbow_y = RElbowAngles_y(idxF2:idxF3);
phase3_RElbow_z = RElbowAngles_z(idxF2:idxF3);

phase4_RElbow_x = RElbowAngles_x(idxF3:idxF4); 
phase4_RElbow_y = RElbowAngles_y(idxF3:idxF4);
phase4_RElbow_z = RElbowAngles_z(idxF3:idxF4);

phase5_RElbow_x = RElbowAngles_x(idxF4:idxF5); 
phase5_RElbow_y = RElbowAngles_y(idxF4:idxF5);
phase5_RElbow_z = RElbowAngles_z(idxF4:idxF5);

% Left Elbow

phase1_LElbow_x = LElbowAngles_x(idxO1:idxF1);
phase1_LElbow_y = LElbowAngles_y(idxO1:idxF1);
phase1_LElbow_z = LElbowAngles_z(idxO1:idxF1);

phase2_LElbow_x = LElbowAngles_x(idxF1:idxF2);
phase2_LElbow_y = LElbowAngles_y(idxF1:idxF2);
phase2_LElbow_z = LElbowAngles_z(idxF1:idxF2);

phase3_LElbow_x = LElbowAngles_x(idxF2:idxF3);
phase3_LElbow_y = LElbowAngles_y(idxF2:idxF3);
phase3_LElbow_z = LElbowAngles_z(idxF2:idxF3);

phase4_LElbow_x = LElbowAngles_x(idxF3:idxF4);
phase4_LElbow_y = LElbowAngles_y(idxF3:idxF4);
phase4_LElbow_z = LElbowAngles_z(idxF3:idxF4);

phase5_LElbow_x = LElbowAngles_x(idxF4:idxF5);
phase5_LElbow_y = LElbowAngles_y(idxF4:idxF5);
phase5_LElbow_z = LElbowAngles_z(idxF4:idxF5);

%% time vectors definition

t2 = 0:1:(length(phase3_RElbow_x)-1);
t3 = 0:1:(length(phase1_RElbow_x)-1);
t4 = 0:1:(length(phase2_RElbow_x)-1);
t5 = 0:1:(length(phase4_RElbow_x)-1);
t6 = 0:1:(length(phase5_RElbow_x)-1);

% Right Elbow figures: MOVEMENT SEGMENTED INTO THE 4 PHASES

figure(10)
subplot(5,1,1)
plot (t3, phase1_RElbow_x, t3, phase1_RElbow_y, '--',t3,phase1_RElbow_z, '.');
legend('x','y','z')
xlabel('Time in frames'); 
ylabel('Degrees');
title('RIGHT ANGLES: Phase 1'); % Picking up the box sagittal plane
subplot(5,1,2)
plot(t4, phase2_RElbow_x,t4,phase2_RElbow_y, '--',t4,phase2_RElbow_z, '.');
legend('x','y','z')
xlabel('Time in frames'); 
ylabel('Degrees');
title('Phase 2'); % Taking steps from sagittal to frontal
subplot(5,1,3)
plot(t2, phase3_RElbow_x, t2, phase3_RElbow_y,'--', t2, phase3_RElbow_z, '.');
legend('x','y','z')
xlabel('Time in frames'); 
ylabel('Degrees');
title('Phase 3');  % Placing the box frontal plane
subplot(5,1,4)
plot(t5, phase4_RElbow_x, t5, phase4_RElbow_y,'--', t5, phase4_RElbow_z, '.');
legend('x','y','z')
xlabel('Time in frames'); 
ylabel('Degrees');
title('Phase 4'); % Taking steps from frontal to sagittal
subplot(5,1,5)
plot(t6, phase5_RElbow_x, t6, phase5_RElbow_y,'--', t6, phase5_RElbow_z, '.');
legend('x','y','z')
xlabel('Time in frames'); 
ylabel('Degrees');
title('Phase 5'); % Placing the box sagittal plane again

% Left Elbow figures: MOVEMENT SEGMENTED INTO THE 4 PHASES

figure(11)
subplot(5,1,1)
plot (t3, phase1_LElbow_x, t3, phase1_LElbow_y, '--',t3,phase1_LElbow_z, '.');
legend('x','y','z')
xlabel('Time in frames'); 
ylabel('Degrees');
title('LEFT ANGLES: Phase 1'); % Picking up the box sagittal plane
subplot(5,1,2)
plot(t4, phase2_LElbow_x,t4, phase2_LElbow_y, '--',t4, phase2_LElbow_z, '.');
legend('x','y','z')
xlabel('Time in frames'); 
ylabel('Degrees');
title('Phase 2'); % Taking steps from sagittal to frontal
subplot(5,1,3)
plot(t2, phase3_LElbow_x, t2, phase3_LElbow_y,'--', t2, phase3_LElbow_z, '.');
legend('x','y','z')
xlabel('Time in frames'); 
ylabel('Degrees');
title('Phase 3'); % Placing the box frontal plane
subplot(5,1,4)
plot(t5, phase4_LElbow_x, t5, phase4_LElbow_y,'--', t5, phase4_LElbow_z, '.');
legend('x','y','z')
xlabel('Time in frames'); 
ylabel('Degrees');
title('Phase 4'); % Taking steps from frontal to sagittal
subplot(5,1,5)
plot(t6, phase5_LElbow_x, t6, phase5_LElbow_y,'--', t6, phase5_LElbow_z, '.');
legend('x','y','z')
xlabel('Time in frames'); 
ylabel('Degrees');
title('Phase 5'); % Placing the box sagittal plane again


%% ROM calculation

% Right Elbow

RElbow_ph1(:,1) = phase1_RElbow_x;
RElbow_ph1(:,2) = phase1_RElbow_y;
RElbow_ph1(:,3) = phase1_RElbow_z;

RElbow_ph2(:,1) = phase2_RElbow_x;
RElbow_ph2(:,2) = phase2_RElbow_y;
RElbow_ph2(:,3) = phase2_RElbow_z;

RElbow_ph3(:,1) = phase3_RElbow_x;
RElbow_ph3(:,2) = phase3_RElbow_y;
RElbow_ph3(:,3) = phase3_RElbow_z;

RElbow_ph4(:,1) = phase4_RElbow_x;
RElbow_ph4(:,2) = phase4_RElbow_y;
RElbow_ph4(:,3) = phase4_RElbow_z;

RElbow_ph5(:,1) = phase5_RElbow_x;
RElbow_ph5(:,2) = phase5_RElbow_y;
RElbow_ph5(:,3) = phase5_RElbow_z;


RElbow_ph1_max = max(RElbow_ph1);
RElbow_ph1_min = min(RElbow_ph1);

RElbow_ph2_max = max(RElbow_ph2);
RElbow_ph2_min = min(RElbow_ph2);

RElbow_ph3_max = max(RElbow_ph3);
RElbow_ph3_min = min(RElbow_ph3);

RElbow_ph4_max = max(RElbow_ph4);
RElbow_ph4_min = min(RElbow_ph4);

RElbow_ph5_max = max(RElbow_ph5);
RElbow_ph5_min = min(RElbow_ph5);

for i = 1:length(RElbow_ph1_max)
    RElbow_ph1_ROM(1,i) = RElbow_ph1_max(i)- RElbow_ph1_min(i);
end

for i = 1:length(RElbow_ph1_max)
    RElbow_ph2_ROM(1,i) = RElbow_ph2_max(i)- RElbow_ph2_min(i);
end

for i = 1:length(RElbow_ph1_max)
    RElbow_ph3_ROM(1,i) = RElbow_ph3_max(i)- RElbow_ph3_min(i);
end

for i = 1:length(RElbow_ph1_max)
    RElbow_ph4_ROM(1,i) = RElbow_ph4_max(i)- RElbow_ph4_min(i);
end

for i = 1:length(RElbow_ph1_max)
    RElbow_ph5_ROM(1,i) = RElbow_ph5_max(i)- RElbow_ph5_min(i);
end

% Left Elbow

LElbow_ph1(:,1) = phase1_LElbow_x;
LElbow_ph1(:,2) = phase1_LElbow_y;
LElbow_ph1(:,3) = phase1_LElbow_z;

LElbow_ph2(:,1) = phase2_LElbow_x;
LElbow_ph2(:,2) = phase2_LElbow_y;
LElbow_ph2(:,3) = phase2_LElbow_z;

LElbow_ph3(:,1) = phase3_LElbow_x;
LElbow_ph3(:,2) = phase3_LElbow_y;
LElbow_ph3(:,3) = phase3_LElbow_z;

LElbow_ph4(:,1) = phase4_LElbow_x;
LElbow_ph4(:,2) = phase4_LElbow_y;
LElbow_ph4(:,3) = phase4_LElbow_z;

LElbow_ph5(:,1) = phase5_LElbow_x;
LElbow_ph5(:,2) = phase5_LElbow_y;
LElbow_ph5(:,3) = phase5_LElbow_z;


LElbow_ph1_max = max(LElbow_ph1);
LElbow_ph1_min = min(LElbow_ph1);

LElbow_ph2_max = max(LElbow_ph2);
LElbow_ph2_min = min(LElbow_ph2);

LElbow_ph3_max = max(LElbow_ph3);
LElbow_ph3_min = min(LElbow_ph3);

LElbow_ph4_max = max(LElbow_ph4);
LElbow_ph4_min = min(LElbow_ph4);

LElbow_ph5_max = max(LElbow_ph5);
LElbow_ph5_min = min(LElbow_ph5);

for i = 1:length(RElbow_ph1_max)
    LElbow_ph1_ROM(1,i) = LElbow_ph1_max(i)- LElbow_ph1_min(i);
end

for i = 1:length(RElbow_ph1_max)
    LElbow_ph2_ROM(1,i) = LElbow_ph2_max(i)- LElbow_ph2_min(i);
end

for i = 1:length(RElbow_ph1_max)
    LElbow_ph3_ROM(1,i) = LElbow_ph3_max(i)- LElbow_ph3_min(i);
end

for i = 1:length(RElbow_ph1_max)
    LElbow_ph4_ROM(1,i) = LElbow_ph4_max(i)- LElbow_ph4_min(i);
end

for i = 1:length(RElbow_ph1_max)
    LElbow_ph5_ROM(1,i) = LElbow_ph5_max(i)- LElbow_ph5_min(i);
end

% Total ROM Elbow Angles

RElbow_total(:,1) = RElbowAngles_x;
RElbow_total(:,2) = RElbowAngles_y;
RElbow_total(:,3) = RElbowAngles_z;

RElbow_total_max = max(RElbow_total);
RElbow_total_min = min(RElbow_total);

for i = 1:length(RElbow_total_max)
    RElbow_total_ROM(1,i) = RElbow_total_max(i)- RElbow_total_min(i);
end

LElbow_total(:,1) = LElbowAngles_x;
LElbow_total(:,2) = LElbowAngles_y;
LElbow_total(:,3) = LElbowAngles_z;

LElbow_total_max = max(LElbow_total);
LElbow_total_min = min(LElbow_total);

for i = 1:length(LElbow_total_max)
    LElbow_total_ROM(1,i) = LElbow_total_max(i)- LElbow_total_min(i);
end

ROM_Right_Elbow = table([RElbow_ph1_ROM(1);RElbow_ph2_ROM(1);RElbow_ph3_ROM(1);RElbow_ph4_ROM(1); RElbow_ph5_ROM(1)],[RElbow_ph1_ROM(2);RElbow_ph2_ROM(2);RElbow_ph3_ROM(2);RElbow_ph4_ROM(2); RElbow_ph5_ROM(2)],[RElbow_ph1_ROM(3);RElbow_ph2_ROM(3);RElbow_ph3_ROM(3);RElbow_ph4_ROM(3); RElbow_ph5_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Phase 1','Phase 2','Phase 3', 'Phase 4', 'Phase 5'})

ROM_Left_Elbow = table([LElbow_ph1_ROM(1);LElbow_ph2_ROM(1);LElbow_ph3_ROM(1);LElbow_ph4_ROM(1); LElbow_ph5_ROM(1)],[LElbow_ph1_ROM(2);LElbow_ph2_ROM(2);LElbow_ph3_ROM(2);LElbow_ph4_ROM(2); LElbow_ph5_ROM(2)],[LElbow_ph1_ROM(3);LElbow_ph2_ROM(3);LElbow_ph3_ROM(3);LElbow_ph4_ROM(3); LElbow_ph5_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Phase 1','Phase 2','Phase 3', 'Phase 4', 'Phase 5'})

TotalROM_Elbow = table([RElbow_total_ROM(1); LElbow_total_ROM(1)],[RElbow_total_ROM(2); LElbow_total_ROM(2)],[RElbow_total_ROM(3); LElbow_total_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Right','Left'})

