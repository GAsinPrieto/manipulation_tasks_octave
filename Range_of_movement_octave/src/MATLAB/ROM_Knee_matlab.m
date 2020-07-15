%% Range of movement of Knee in Sagittal Lifting

% Code AUTHOR: Yaiza Benito Molpeceres. DATE: January-May 2020.

% Segmentation is based on box signal and consists of five phases shown in figure 1:
    % 1: Subject starts at upright position and lowers down to pick the box (idxO1:idxF1)
    % 2: Subject holds and deposits the box (idxF1:idxF2)
    % 3: Subject raises up to upright position and prepares to take the box again (idxF2:idxO3)
    % 4: Subject holds and deposits the box in its original location (idxO3:idxF3)
    % 5: Subject raises up to upright position (idxF3:idxF4)  
    
% ROM is calculated for phases 1, 2, 4 and 5, which are renamed as 1, 2, 3 and 4.clear all % Clear variables

clear all % Clear variables
close all % Close figures
clc

load('..\tests\data\input\dinamica59_B.mat')
             
Ts = 1/Fs;
t_total = (double(frames)*Ts);
%t_1000 = 0:(t_total/1000):(t_total-(t_total/1000));
t_100 = 1:1:100;

% Normal Ranges of Knee Motion
for i =1:(length(t_100))
    flex_knee(1,i) = 130; % flexion
    ext_knee(1,i) = -15;   % extension
    intRot_hip(1,i) = 10;   % internal rotation
end

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
%idxF1 = idx_03(1)-1;                               % solo para dinamica 03
idxF2 = idx(idx_flat2(1)); 
idxO3 = idx(idx_flat3(length(idx_flat3))); 
idxF3 = idx(idx_flat4(3));                         % para dinamica 59 es 4, para din 03 es 2 y para el resto es 1                    
idxF4 = length(position_Z);
 
% Visualization of the segmentation
figure(1)
plot(t,position_Z,t(locs),position_Z(locs),'o');
hold on
plot(t,position_Z,'r*','MarkerIndices',idxO1);
plot(t,position_Z,'r*','MarkerIndices',idxF1);
plot(t,position_Z,'r*','MarkerIndices',idxF2);
plot(t,position_Z,'r*','MarkerIndices',idxO3);
plot(t,position_Z,'r*','MarkerIndices',idxF3);
plot(t,position_Z,'r*','MarkerIndices',idxF4);
hold off
xlabel('Time in seconds');
ylabel('Millimetres');
title('Box signal segmentation in 4 phases'); % or title: Sagittal Lifting Segmentation

%% Import model data for Knee Relative Angles

        RKneeAngles_x = ModelData.Raw.(ModelOutput{3})(1,:); 
        RKneeAngles_y = ModelData.Raw.(ModelOutput{3})(2,:);
        RKneeAngles_z = ModelData.Raw.(ModelOutput{3})(3,:);
       
        LKneeAngles_x = ModelData.Raw.(ModelOutput{4})(1,:);
        LKneeAngles_y = ModelData.Raw.(ModelOutput{4})(2,:);
        LKneeAngles_z = ModelData.Raw.(ModelOutput{4})(3,:);
       
        
%% Knee signals visualization before going through segmentation

figure(2)
frames = length(position_Z);
% (1:(frames/100):frames)  si frames>1000
% (1:(frames/100):frames+(frames/100)) si frames<1000

plot(t_100,RKneeAngles_x(1:(frames/100):frames)); 
hold on 
plot(t_100, flex_knee,'--',t_100, ext_knee,'--'); 
legend({'x', 'Normal range of flexion-extension'},'Location','best');
xlabel('% Lifting Cycle'); 
ylabel('Degrees');
%ylim([-100 160]);
title('Right Knee angles: Flexo-Extension');

figure(3)
plot(t_100, RKneeAngles_y(1:(frames/100):frames));
legend({'y'},'Location','best');
xlabel('% Lifting Cycle'); 
ylabel('Degrees');  
%ylim([-70 50]);
title('Right Knee angles: Abduction-Adduction');

figure(4)
plot(t_100, RKneeAngles_z(1:(frames/100):frames));
hold on 
plot(t_100, intRot_hip,'--'); 
legend({'z', 'Normal range of internal rotation'},'Location','best');
xlabel('% Lifting Cycle');   
ylabel('Degrees');
ylim([-70 20]);
title('Right Knee angles: Internal-External Rotation');

figure(5)
% (1:(frames/100):frames)  si frames>1000
% (1:(frames/100):frames+(frames/100)) si frames<1000

plot(t_100,LKneeAngles_x(1:(frames/100):frames)); 
hold on 
plot(t_100, flex_knee,'--',t_100, ext_knee,'--'); 
legend({'x', 'Normal range of flexion-extension'},'Location','best');
xlabel('% Lifting Cycle');  
ylabel('Degrees');
%ylim([-100 160]);
title('Left Knee angles: Flexo-Extension');

figure(6)
plot(t_100, LKneeAngles_y(1:(frames/100):frames));
legend({'y'},'Location','best');
xlabel('% Lifting Cycle'); 
ylabel('Degrees');
%ylim([-70 50]);
title('Left Knee angles: Abduction-Adduction');

figure(7)
plot(t_100, LKneeAngles_z(1:(frames/100):frames));
hold on 
plot(t_100, intRot_hip,'--'); 
legend({'z', 'Normal range of internal rotation'},'Location','best');
xlabel('% Lifting Cycle'); 
ylabel('Degrees');
%ylim([-70 80]);
title('Left Knee angles: Internal-External Rotation');

%figure(3)
%plot(t ,RKneeAngles_x, t, RKneeAngles_y, t, RKneeAngles_z);
%legend('x','y','z');
%xlabel('time in frames'); % frames totales del trial
%title('Right Knee angles');

%figure(4)
%plot(t,LKneeAngles_x, t, LKneeAngles_y, t, LKneeAngles_z);   
%legend('x','y','z');
%xlabel('time in frames'); % frames totales del trial
%title('Left Knee angles');

%% Signal Segmentation

% Right Knee

phase1_RKnee_x = RKneeAngles_x(idxO1:idxF1);  
phase1_RKnee_y = RKneeAngles_y(idxO1:idxF1);
phase1_RKnee_z = RKneeAngles_z(idxO1:idxF1);

phase2_RKnee_x = RKneeAngles_x(idxF1:idxF2); 
phase2_RKnee_y = RKneeAngles_y(idxF1:idxF2);
phase2_RKnee_z = RKneeAngles_z(idxF1:idxF2);

phase3_RKnee_x = RKneeAngles_x(idxO3:idxF3); 
phase3_RKnee_y = RKneeAngles_y(idxO3:idxF3);
phase3_RKnee_z = RKneeAngles_z(idxO3:idxF3);

phase4_RKnee_x = RKneeAngles_x(idxF3:idxF4); 
phase4_RKnee_y = RKneeAngles_y(idxF3:idxF4);
phase4_RKnee_z = RKneeAngles_z(idxF3:idxF4);

% Left Knee

phase1_LKnee_x = LKneeAngles_x(idxO1:idxF1);
phase1_LKnee_y = LKneeAngles_y(idxO1:idxF1);
phase1_LKnee_z = LKneeAngles_z(idxO1:idxF1);

phase2_LKnee_x = LKneeAngles_x(idxF1:idxF2);
phase2_LKnee_y = LKneeAngles_y(idxF1:idxF2);
phase2_LKnee_z = LKneeAngles_z(idxF1:idxF2);

phase3_LKnee_x = LKneeAngles_x(idxO3:idxF3);
phase3_LKnee_y = LKneeAngles_y(idxO3:idxF3);
phase3_LKnee_z = LKneeAngles_z(idxO3:idxF3);

phase4_LKnee_x = LKneeAngles_x(idxF3:idxF4);
phase4_LKnee_y = LKneeAngles_y(idxF3:idxF4);
phase4_LKnee_z = LKneeAngles_z(idxF3:idxF4);

%% time vectors definition

t2 = (0:1:(length(phase3_RKnee_x)-1))*Ts;
t3 = (0:1:(length(phase1_RKnee_x)-1))*Ts;
t4 = (0:1:(length(phase2_RKnee_x)-1))*Ts;
t5 = (0:1:(length(phase4_RKnee_x)-1))*Ts;

% Right Knee figures: MOVEMENT SEGMENTED INTO THE 4 PHASES

figure(8)
subplot(4,1,1)
plot (t3, phase1_RKnee_x, t3, phase1_RKnee_y, '--',t3,phase1_RKnee_z, '.');
legend('x','y','z')
xlabel('Time in seconds');
ylabel('Degrees');
title('RIGHT ANGLES: Phase 1');  % Lowering without box Right Knee
subplot(4,1,2)
plot(t4, phase2_RKnee_x,t4,phase2_RKnee_y, '--',t4,phase2_RKnee_z, '.');
legend('x','y','z')
xlabel('Time in seconds');
ylabel('Degrees');
title('Phase 2');  % Raising up with box Right Knee
subplot(4,1,3)
plot(t2, phase3_RKnee_x, t2, phase3_RKnee_y,'--', t2, phase3_RKnee_z, '.');
legend('x','y','z')
xlabel('Time in seconds');
ylabel('Degrees');
title('Phase 3');  % Lowering with box Right Knee
subplot(4,1,4)
plot(t5, phase4_RKnee_x, t5, phase4_RKnee_y,'--', t5, phase4_RKnee_z, '.');
legend('x','y','z')
xlabel('Time in seconds');
ylabel('Degrees');
title('Phase 4');  % Raising up without box Right Knee

% Left Knee figures: MOVEMENT SEGMENTED INTO THE 4 PHASES

figure(9)
subplot(4,1,1)
plot (t3, phase1_LKnee_x, t3, phase1_LKnee_y, '--',t3,phase1_LKnee_z, '.');
legend('x','y','z')
xlabel('Time in seconds');
ylabel('Degrees');
title('LEFT ANGLES: Phase 1');  % Lowering without box Right Knee
subplot(4,1,2)
plot(t4, phase2_LKnee_x,t4, phase2_LKnee_y, '--',t4, phase2_LKnee_z, '.');
legend('x','y','z')
xlabel('Time in seconds');
ylabel('Degrees');
title('Phase 2');  % Raising up with box Right Knee
subplot(4,1,3)
plot(t2, phase3_LKnee_x, t2, phase3_LKnee_y,'--', t2, phase3_LKnee_z, '.');
legend('x','y','z')
xlabel('Time in seconds');
ylabel('Degrees');
title('Phase 3');  % Lowering with box Right Knee
subplot(4,1,4)
plot(t5, phase4_LKnee_x, t5, phase4_LKnee_y,'--', t5, phase4_LKnee_z, '.');
legend('x','y','z')
xlabel('Time in seconds');
ylabel('Degrees');
title('Phase 4');  % Raising up without box Right Knee

%% ROM calculation

% Right Knee

RKnee_ph1(:,1) = phase1_RKnee_x;
RKnee_ph1(:,2) = phase1_RKnee_y;
RKnee_ph1(:,3) = phase1_RKnee_z;

RKnee_ph2(:,1) = phase2_RKnee_x;
RKnee_ph2(:,2) = phase2_RKnee_y;
RKnee_ph2(:,3) = phase2_RKnee_z;

RKnee_ph3(:,1) = phase3_RKnee_x;
RKnee_ph3(:,2) = phase3_RKnee_y;
RKnee_ph3(:,3) = phase3_RKnee_z;

RKnee_ph4(:,1) = phase4_RKnee_x;
RKnee_ph4(:,2) = phase4_RKnee_y;
RKnee_ph4(:,3) = phase4_RKnee_z;

RKnee_ph1_max = max(RKnee_ph1);
RKnee_ph1_min = min(RKnee_ph1);

RKnee_ph2_max = max(RKnee_ph2);
RKnee_ph2_min = min(RKnee_ph2);

RKnee_ph3_max = max(RKnee_ph3);
RKnee_ph3_min = min(RKnee_ph3);

RKnee_ph4_max = max(RKnee_ph4);
RKnee_ph4_min = min(RKnee_ph4);

for i = 1:length(RKnee_ph1_max)
    RKnee_ph1_ROM(1,i) = RKnee_ph1_max(i)- RKnee_ph1_min(i);
end

for i = 1:length(RKnee_ph1_max)
    RKnee_ph2_ROM(1,i) = RKnee_ph2_max(i)- RKnee_ph2_min(i);
end

for i = 1:length(RKnee_ph1_max)
    RKnee_ph3_ROM(1,i) = RKnee_ph3_max(i)- RKnee_ph3_min(i);
end

for i = 1:length(RKnee_ph1_max)
    RKnee_ph4_ROM(1,i) = RKnee_ph4_max(i)- RKnee_ph4_min(i);
end

% Left Knee

LKnee_ph1(:,1) = phase1_LKnee_x;
LKnee_ph1(:,2) = phase1_LKnee_y;
LKnee_ph1(:,3) = phase1_LKnee_z;

LKnee_ph2(:,1) = phase2_LKnee_x;
LKnee_ph2(:,2) = phase2_LKnee_y;
LKnee_ph2(:,3) = phase2_LKnee_z;

LKnee_ph3(:,1) = phase3_LKnee_x;
LKnee_ph3(:,2) = phase3_LKnee_y;
LKnee_ph3(:,3) = phase3_LKnee_z;

LKnee_ph4(:,1) = phase4_LKnee_x;
LKnee_ph4(:,2) = phase4_LKnee_y;
LKnee_ph4(:,3) = phase4_LKnee_z;

LKnee_ph1_max = max(LKnee_ph1);
LKnee_ph1_min = min(LKnee_ph1);

LKnee_ph2_max = max(LKnee_ph2);
LKnee_ph2_min = min(LKnee_ph2);

LKnee_ph3_max = max(LKnee_ph3);
LKnee_ph3_min = min(LKnee_ph3);

LKnee_ph4_max = max(LKnee_ph4);
LKnee_ph4_min = min(LKnee_ph4);

for i = 1:length(LKnee_ph1_max)
    LKnee_ph1_ROM(1,i) = LKnee_ph1_max(i)- LKnee_ph1_min(i);
end

for i = 1:length(LKnee_ph1_max)
    LKnee_ph2_ROM(1,i) = LKnee_ph2_max(i)- LKnee_ph2_min(i);
end

for i = 1:length(LKnee_ph1_max)
    LKnee_ph3_ROM(1,i) = LKnee_ph3_max(i)- LKnee_ph3_min(i);
end

for i = 1:length(LKnee_ph1_max)
    LKnee_ph4_ROM(1,i) = LKnee_ph4_max(i)- LKnee_ph4_min(i);
end

% Total ROM Knee Angles

RKnee_total(:,1) = RKneeAngles_x;
RKnee_total(:,2) = RKneeAngles_y;
RKnee_total(:,3) = RKneeAngles_z;

RKnee_total_max = max(RKnee_total);
RKnee_total_min = min(RKnee_total);

for i = 1:length(RKnee_total_max)
    RKnee_total_ROM(1,i) = RKnee_total_max(i)- RKnee_total_min(i);
end

LKnee_total(:,1) = LKneeAngles_x;
LKnee_total(:,2) = LKneeAngles_y;
LKnee_total(:,3) = LKneeAngles_z;

LKnee_total_max = max(LKnee_total);
LKnee_total_min = min(LKnee_total);

for i = 1:length(LKnee_total_max)
    LKnee_total_ROM(1,i) = LKnee_total_max(i)- LKnee_total_min(i);
end

ROM_Right_Knee = table([RKnee_ph1_ROM(1);RKnee_ph2_ROM(1);RKnee_ph3_ROM(1);RKnee_ph4_ROM(1)],[RKnee_ph1_ROM(2);RKnee_ph2_ROM(2);RKnee_ph3_ROM(2);RKnee_ph4_ROM(2)],[RKnee_ph1_ROM(3);RKnee_ph2_ROM(3);RKnee_ph3_ROM(3);RKnee_ph4_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Phase 1','Phase 2','Phase 3', 'Phase 4'})

ROM_Left_Knee = table([LKnee_ph1_ROM(1);LKnee_ph2_ROM(1);LKnee_ph3_ROM(1);LKnee_ph4_ROM(1)],[LKnee_ph1_ROM(2);LKnee_ph2_ROM(2);LKnee_ph3_ROM(2);LKnee_ph4_ROM(2)],[LKnee_ph1_ROM(3);LKnee_ph2_ROM(3);LKnee_ph3_ROM(3);LKnee_ph4_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Phase 1','Phase 2','Phase 3', 'Phase 4'})

TotalROM_Knee = table([RKnee_total_ROM(1); LKnee_total_ROM(1)],[RKnee_total_ROM(2); LKnee_total_ROM(2)],[RKnee_total_ROM(3); LKnee_total_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Right','Left'})