%% Range of movement of Spine in Sagittal Lifting

% Code AUTHOR: Yaiza Benito Molpeceres. DATE: January-May 2020.
% Adapted to Octave by Guillermo Asín Prieto

% Segmentation is based on box signal and consists of five phases shown in figure 1:
    % 1: Subject starts at upright position and lowers down to pick the box (idxO1:idxF1)
    % 2: Subject holds and deposits the box (idxF1:idxF2)
    % 3: Subject raises up to upright position and prepares to take the box again (idxF2:idxO3)
    % 4: Subject holds and deposits the box in its original location (idxO3:idxF3)
    % 5: Subject raises up to upright position (idxF3:idxF4)  
    
% ROM is calculated for phases 1, 2, 4 and 5, which are renamed as 1, 2, 3 and 4.

clear all % Clear variables
close all % Close figures
clc

load('dinamica44_B.mat')
          
Ts = 1/Fs;
t_total = (double(frames)*Ts);
%t_1000 = 0:(t_total/1000):(t_total-(t_total/1000));
t_100 = 1:1:100;

% Normal Ranges of Spine Motion
for i =1:(length(t_100))
    flex_spine(1,i) = 75; % flexion
    ext_spine(1,i) = -30;   % extension
    abd_spine(1,i) = -35;   % abduction
end

%% Events of interest in vertical box trajectory identification
 
position_Z = SIDEBOX3(:,3)';
 
[pks, locs] = findpeaks(position_Z,'minPeakProminence',20);
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
idxF3 = idx(idx_flat4(1));                         % para dinamica 59 es 4, para din 03 y din 61 es 2 y para el resto es 1                    
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

%% Import model data for Spine Relative Angles

        RSpineAngles_x = ModelData.Raw.(ModelOutput{11})(1,:); 
        RSpineAngles_y = ModelData.Raw.(ModelOutput{11})(2,:);
        RSpineAngles_z = ModelData.Raw.(ModelOutput{11})(3,:);
        

        LSpineAngles_x = ModelData.Raw.(ModelOutput{12})(1,:);
        LSpineAngles_y = ModelData.Raw.(ModelOutput{12})(2,:);
        LSpineAngles_z = ModelData.Raw.(ModelOutput{12})(3,:);
        
        
%% Spine signals visualization before going through segmentation

figure(2)
frames = length(position_Z);
% (1:(frames/100):frames)  si frames>1000
% (1:(frames/100):frames+(frames/100)) si frames<1000

plot(t_100,RSpineAngles_x(1:(frames/100):frames)); 
hold on 
plot(t_100, flex_spine,'--',t_100, ext_spine,'--'); 
legend({'x', 'Normal range of flexion-extension'},'Location','best');
xlabel('% Lifting Cycle');
ylabel('Degrees');
%ylim([-100 160]);
title('Right Spine angles: Flexo-Extension');

figure(3)
plot(t_100, RSpineAngles_y(1:(frames/100):frames));
hold on 
plot(t_100, abd_spine,'--'); 
legend({'y', 'Normal range of abduction'},'Location','best');
xlabel('% Lifting Cycle');
ylabel('Degrees');
ylim([-40 10]);
title('Right Spine angles: Abduction-Adduction');

figure(4)
plot(t_100, RSpineAngles_z(1:(frames/100):frames));
legend({'z'},'Location','best');
xlabel('% Lifting Cycle');
ylabel('Degrees');
%ylim([-70 50]);
title('Right Spine angles: Internal-External Rotation');

figure(5)
% (1:(frames/100):frames)  si frames>1000
% (1:(frames/100):frames+(frames/100)) si frames<1000

plot(t_100,LSpineAngles_x(1:(frames/100):frames)); 
hold on 
plot(t_100, flex_spine,'--',t_100, ext_spine,'--'); 
legend({'x', 'Normal range of flexion-extension'},'Location','best');
xlabel('% Lifting Cycle');
ylabel('Degrees');
%ylim([-100 160]);
title('Left Spine angles: Flexo-Extension');

figure(6)
plot(t_100, LSpineAngles_y(1:(frames/100):frames));
hold on 
plot(t_100, abd_spine,'--'); 
legend({'y', 'Normal range of abduction'},'Location','best');
xlabel('% Lifting Cycle');
ylabel('Degrees');
ylim([-40 5]);
title('Left Spine angles: Abduction-Adduction');

figure(7)
plot(t_100, LSpineAngles_z(1:(frames/100):frames));
legend({'z'},'Location','best');
xlabel('% Lifting Cycle');
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

% Right Spine (I keep distinguishing right and left because although it makes no difference for x axis, it does for y and z ones). 

phase1_RSpine_x = RSpineAngles_x(idxO1:idxF1);  
phase1_RSpine_y = RSpineAngles_y(idxO1:idxF1);
phase1_RSpine_z = RSpineAngles_z(idxO1:idxF1);

phase2_RSpine_x = RSpineAngles_x(idxF1:idxF2); 
phase2_RSpine_y = RSpineAngles_y(idxF1:idxF2);
phase2_RSpine_z = RSpineAngles_z(idxF1:idxF2);

phase3_RSpine_x = RSpineAngles_x(idxO3:idxF3); 
phase3_RSpine_y = RSpineAngles_y(idxO3:idxF3);
phase3_RSpine_z = RSpineAngles_z(idxO3:idxF3);

phase4_RSpine_x = RSpineAngles_x(idxF3:idxF4); 
phase4_RSpine_y = RSpineAngles_y(idxF3:idxF4);
phase4_RSpine_z = RSpineAngles_z(idxF3:idxF4);

% Left Spine

phase1_LSpine_x = LSpineAngles_x(idxO1:idxF1);
phase1_LSpine_y = LSpineAngles_y(idxO1:idxF1);
phase1_LSpine_z = LSpineAngles_z(idxO1:idxF1);

phase2_LSpine_x = LSpineAngles_x(idxF1:idxF2);
phase2_LSpine_y = LSpineAngles_y(idxF1:idxF2);
phase2_LSpine_z = LSpineAngles_z(idxF1:idxF2);

phase3_LSpine_x = LSpineAngles_x(idxO3:idxF3);
phase3_LSpine_y = LSpineAngles_y(idxO3:idxF3);
phase3_LSpine_z = LSpineAngles_z(idxO3:idxF3);

phase4_LSpine_x = LSpineAngles_x(idxF3:idxF4);
phase4_LSpine_y = LSpineAngles_y(idxF3:idxF4);
phase4_LSpine_z = LSpineAngles_z(idxF3:idxF4);

%% time vectors definition

t2 = (0:1:(length(phase3_RSpine_x)-1))*Ts;
t3 = (0:1:(length(phase1_RSpine_x)-1))*Ts;
t4 = (0:1:(length(phase2_RSpine_x)-1))*Ts;
t5 = (0:1:(length(phase4_RSpine_x)-1))*Ts;

% Right Spine figures: MOVEMENT SEGMENTED INTO THE 4 PHASES

figure(8)
subplot(4,1,1)
plot (t3, phase1_RSpine_x, t3, phase1_RSpine_y, '--',t3,phase1_RSpine_z, '.');
legend('x','y','z')
xlabel('Time in seconds'); 
ylabel('Degrees');
title('RIGHT ANGLES: Phase 1');  % Lowering without box Right Spine
subplot(4,1,2)
plot(t4, phase2_RSpine_x,t4,phase2_RSpine_y, '--',t4,phase2_RSpine_z, '.');
legend('x','y','z')
xlabel('Time in seconds');  
ylabel('Degrees');
title('Phase 2');  %Raising up with box Right Spine
subplot(4,1,3)
plot(t2, phase3_RSpine_x, t2, phase3_RSpine_y,'--', t2, phase3_RSpine_z, '.');
legend('x','y','z')
xlabel('Time in seconds');
ylabel('Degrees');
title('Phase 3');  % Lowering with box Right Spine
subplot(4,1,4)
plot(t5, phase4_RSpine_x, t5, phase4_RSpine_y,'--', t5, phase4_RSpine_z, '.');
xlabel('Time in seconds');  
ylabel('Degrees');
title('Phase 4');  % Raising up without box Right Spine

% Left Spine figures: MOVEMENT SEGMENTED INTO THE 4 PHASES

figure(9)
subplot(4,1,1)
plot (t3, phase1_LSpine_x, t3, phase1_LSpine_y, '--',t3,phase1_LSpine_z, '.');
legend('x','y','z')
xlabel('Time in seconds');  
ylabel('Degrees');
title('LEFT ANGLES: Phase 1');  % Lowering without box Right Spine
subplot(4,1,2)
plot(t4, phase2_LSpine_x,t4, phase2_LSpine_y, '--',t4, phase2_LSpine_z, '.');
legend('x','y','z')
xlabel('Time in seconds'); 
ylabel('Degrees');
title('Phase 2');  %Raising up with box Left Spine
subplot(4,1,3)
plot(t2, phase3_LSpine_x, t2, phase3_LSpine_y,'--', t2, phase3_LSpine_z, '.');
legend('x','y','z')
xlabel('Time in seconds');  
ylabel('Degrees');
title('Phase 3');  % Lowering with box Left Spine
subplot(4,1,4)
plot(t5, phase4_LSpine_x, t5, phase4_LSpine_y,'--', t5, phase4_LSpine_z, '.');
legend('x','y','z')
xlabel('Time in seconds'); 
ylabel('Degrees');
title('Phase 4');  % Raising up without box Left Spine

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

RSpine_ph4(:,1) = phase4_RSpine_x;
RSpine_ph4(:,2) = phase4_RSpine_y;
RSpine_ph4(:,3) = phase4_RSpine_z;

RSpine_ph1_max = max(RSpine_ph1);
RSpine_ph1_min = min(RSpine_ph1);

RSpine_ph2_max = max(RSpine_ph2);
RSpine_ph2_min = min(RSpine_ph2);

RSpine_ph3_max = max(RSpine_ph3);
RSpine_ph3_min = min(RSpine_ph3);

RSpine_ph4_max = max(RSpine_ph4);
RSpine_ph4_min = min(RSpine_ph4);

for i = 1:length(RSpine_ph1_max)
    RSpine_ph1_ROM(1,i) = RSpine_ph1_max(i)- RSpine_ph1_min(i);
end

for i = 1:length(RSpine_ph1_max)
    RSpine_ph2_ROM(1,i) = RSpine_ph2_max(i)- RSpine_ph2_min(i);
end

for i = 1:length(RSpine_ph1_max)
    RSpine_ph3_ROM(1,i) = RSpine_ph3_max(i)- RSpine_ph3_min(i);
end

for i = 1:length(RSpine_ph1_max)
    RSpine_ph4_ROM(1,i) = RSpine_ph4_max(i)- RSpine_ph4_min(i);
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

LSpine_ph4(:,1) = phase4_LSpine_x;
LSpine_ph4(:,2) = phase4_LSpine_y;
LSpine_ph4(:,3) = phase4_LSpine_z;

LSpine_ph1_max = max(LSpine_ph1);
LSpine_ph1_min = min(LSpine_ph1);

LSpine_ph2_max = max(LSpine_ph2);
LSpine_ph2_min = min(LSpine_ph2);

LSpine_ph3_max = max(LSpine_ph3);
LSpine_ph3_min = min(LSpine_ph3);

LSpine_ph4_max = max(LSpine_ph4);
LSpine_ph4_min = min(LSpine_ph4);

for i = 1:length(RSpine_ph1_max)
    LSpine_ph1_ROM(1,i) = LSpine_ph1_max(i)- LSpine_ph1_min(i);
end

for i = 1:length(RSpine_ph1_max)
    LSpine_ph2_ROM(1,i) = LSpine_ph2_max(i)- LSpine_ph2_min(i);
end

for i = 1:length(RSpine_ph1_max)
    LSpine_ph3_ROM(1,i) = LSpine_ph3_max(i)- LSpine_ph3_min(i);
end

for i = 1:length(RSpine_ph1_max)
    LSpine_ph4_ROM(1,i) = LSpine_ph4_max(i)- LSpine_ph4_min(i);
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

ROM_Right_Spine = table([RSpine_ph1_ROM(1);RSpine_ph2_ROM(1);RSpine_ph3_ROM(1);RSpine_ph4_ROM(1)],[RSpine_ph1_ROM(2);RSpine_ph2_ROM(2);RSpine_ph3_ROM(2);RSpine_ph4_ROM(2)],[RSpine_ph1_ROM(3);RSpine_ph2_ROM(3);RSpine_ph3_ROM(3);RSpine_ph4_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Phase 1','Phase 2','Phase 3', 'Phase 4'})

ROM_Left_Spine = table([LSpine_ph1_ROM(1);LSpine_ph2_ROM(1);LSpine_ph3_ROM(1);LSpine_ph4_ROM(1)],[LSpine_ph1_ROM(2);LSpine_ph2_ROM(2);LSpine_ph3_ROM(2);LSpine_ph4_ROM(2)],[LSpine_ph1_ROM(3);LSpine_ph2_ROM(3);LSpine_ph3_ROM(3);LSpine_ph4_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Phase 1','Phase 2','Phase 3', 'Phase 4'})

TotalROM_Spine = table([RSpine_total_ROM(1); LSpine_total_ROM(1)],[RSpine_total_ROM(2); LSpine_total_ROM(2)],[RSpine_total_ROM(3); LSpine_total_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Right','Left'})
