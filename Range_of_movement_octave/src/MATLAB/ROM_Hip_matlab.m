%% Range of movement of Hip in Sagittal Lifting

% Code AUTHOR: Yaiza Benito Molpeceres. DATE: January-May 2020.

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

load('..\tests\data\input\Dinamica61_B.mat')
           
Ts = 1/Fs;                              
t_total = (double(frames)*Ts);
%t_1000 = 0:(t_total/1000):(t_total-(t_total/1000));
t_100 = 1:1:100;

% Normal Ranges of Hip Motion
for i =1:(length(t_100))
    flex_hip(1,i) = 130; % flexion
    ext_hip(1,i) = -30;   % extension
    abd_hip(1,i) = -50;   % abduction
    add_hip(1,i) = 30;   % adduction
    intRot_hip(1,i) = 40;   % internal rotation
    extRot_hip(1,i) = -45;   % external rotation
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
idxF3 = idx(idx_flat4(1));                         % para dinamica 59 es 3, para din 03 es 2 y para el resto es 1                    
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
title('Box signal segmentation in 4 phases');

%% Import model data for Hip Relative Angles

        RHipAngles_x = ModelData.Raw.(ModelOutput{1})(1,:); 
        RHipAngles_y = ModelData.Raw.(ModelOutput{1})(2,:);
        RHipAngles_z = ModelData.Raw.(ModelOutput{1})(3,:);
       
        LHipAngles_x = ModelData.Raw.(ModelOutput{2})(1,:);
        LHipAngles_y = ModelData.Raw.(ModelOutput{2})(2,:);
        LHipAngles_z = ModelData.Raw.(ModelOutput{2})(3,:);
       

%% Hip signals visualization before going through segmentation

figure(2)
frames = length(position_Z);
% (1:(frames/100):frames)  si frames>1000
% (1:(frames/100):frames+(frames/100)) si frames<1000

plot(t_100,RHipAngles_x(1:(frames/100):frames)); 
hold on 
plot(t_100, flex_hip,'--',t_100, ext_hip,'--'); 
legend({'x', 'Normal range of flexion-extension'},'Location','best');
xlabel('% Lifting Cycle');
ylabel('Degrees');
ylim([-100 160]);
title('Right Hip angles: Flexo-Extension');

figure(3)
plot(t_100, RHipAngles_y(1:(frames/100):frames));
hold on 
plot(t_100, abd_hip,'--',t_100, add_hip,'--'); 
legend({'y', 'Normal range of abduction-adduction'},'Location','best');
xlabel('% Lifting Cycle)');  
ylabel('Degrees');
ylim([-70 50]);
title('Right Hip angles: Abduction-Adduction');

figure(4)
plot(t_100, RHipAngles_z(1:(frames/100):frames));
hold on 
plot(t_100, intRot_hip,'--',t_100, extRot_hip,'--'); 
legend({'z', 'Normal range of internal-external rotation'},'Location','best');
xlabel('% Lifting Cycle');  
ylabel('Degrees');
ylim([-70 50]);
title('Right Hip angles: Internal-External Rotation');

figure(5)
% (1:(frames/100):frames)  si frames>1000
% (1:(frames/100):frames+(frames/100)) si frames<1000

plot(t_100,LHipAngles_x(1:(frames/100):frames)); 
hold on 
plot(t_100, flex_hip,'--',t_100, ext_hip,'--'); 
legend({'x', 'Normal range of flexion-extension'},'Location','best');
xlabel('% Lifting Cycle');
ylabel('Degrees');
ylim([-100 160]);
title('Left Hip angles: Flexo-Extension');

figure(6)
plot(t_100, LHipAngles_y(1:(frames/100):frames));
hold on 
plot(t_100, abd_hip,'--',t_100, add_hip,'--'); 
legend({'y', 'Normal range of abduction-adduction'},'Location','best');
xlabel('% Lifting Cycle'); 
ylabel('Degrees');
ylim([-70 50]);
title('Left Hip angles: Abduction-Adduction');

figure(7)
plot(t_100, LHipAngles_z(1:(frames/100):frames));
hold on 
plot(t_100, intRot_hip,'--',t_100, extRot_hip,'--'); 
legend({'z', 'Normal range of internal-external rotation'},'Location','best');
xlabel('% Lifting Cycle');
ylabel('Degrees');
ylim([-70 80]);
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

phase2_RHip_x = RHipAngles_x(idxF1:idxF2); 
phase2_RHip_y = RHipAngles_y(idxF1:idxF2);
phase2_RHip_z = RHipAngles_z(idxF1:idxF2);

phase3_RHip_x = RHipAngles_x(idxO3:idxF3); 
phase3_RHip_y = RHipAngles_y(idxO3:idxF3);
phase3_RHip_z = RHipAngles_z(idxO3:idxF3);

phase4_RHip_x = RHipAngles_x(idxF3:idxF4); 
phase4_RHip_y = RHipAngles_y(idxF3:idxF4);
phase4_RHip_z = RHipAngles_z(idxF3:idxF4);

% Left Hip 

phase1_LHip_x = LHipAngles_x(idxO1:idxF1);
phase1_LHip_y = LHipAngles_y(idxO1:idxF1);
phase1_LHip_z = LHipAngles_z(idxO1:idxF1);

phase2_LHip_x = LHipAngles_x(idxF1:idxF2);
phase2_LHip_y = LHipAngles_y(idxF1:idxF2);
phase2_LHip_z = LHipAngles_z(idxF1:idxF2);

phase3_LHip_x = LHipAngles_x(idxO3:idxF3);
phase3_LHip_y = LHipAngles_y(idxO3:idxF3);
phase3_LHip_z = LHipAngles_z(idxO3:idxF3);

phase4_LHip_x = LHipAngles_x(idxF3:idxF4);
phase4_LHip_y = LHipAngles_y(idxF3:idxF4);
phase4_LHip_z = LHipAngles_z(idxF3:idxF4);

%% time vectors definition

t2 = (0:1:length(phase3_RHip_x)-1)*Ts;
t3 = (0:1:(length(phase1_RHip_x)-1))*Ts;
t4 = (0:1:(length(phase2_RHip_x)-1))*Ts;
t5 = (0:1:(length(phase4_RHip_x)-1))*Ts;

% Right Hip figures: MOVEMENT SEGMENTED INTO THE 4 PHASES

figure(8)
subplot(4,1,1)
plot (t3, phase1_RHip_x, t3, phase1_RHip_y, '--',t3,phase1_RHip_z, '.');
legend('x','y','z')
xlabel('Time in seconds'); 
ylabel('Degrees');
title('RIGHT ANGLES: Phase 1');  % Lowering without box Right Hip
subplot(4,1,2)
plot(t4, phase2_RHip_x,t4,phase2_RHip_y, '--',t4,phase2_RHip_z, '.');
legend('x','y','z')
xlabel('Time in seconds'); 
ylabel('Degrees');
title('Phase 2');  %Raising up with box Right Hip
subplot(4,1,3)
plot(t2, phase3_RHip_x, t2, phase3_RHip_y,'--', t2, phase3_RHip_z, '.');
legend('x','y','z')
xlabel('Time in seconds');  
ylabel('Degrees');
title('Phase 3');  % Lowering with box Right Hip
subplot(4,1,4)
plot(t5, phase4_RHip_x, t5, phase4_RHip_y,'--', t5, phase4_RHip_z, '.');
legend('x','y','z')
xlabel('Time in seconds');  
ylabel('Degrees');
title('Phase 4');  % Raising up without box Right Hip

% Left Hip figures: MOVEMENT SEGMENTED INTO THE 4 PHASES

figure(9)
subplot(4,1,1)
plot (t3, phase1_LHip_x, t3, phase1_LHip_y, '--',t3,phase1_LHip_z, '.');
legend('x','y','z')
xlabel('Time in seconds'); 
ylabel('Degrees');
title('LEFT ANGLES: Phase 1');  % Lowering without box Right Hip
subplot(4,1,2)
plot(t4, phase2_LHip_x,t4, phase2_LHip_y, '--',t4, phase2_LHip_z, '.');
legend('x','y','z')
xlabel('Time in seconds'); 
ylabel('Degrees');
title('Phase 2');  %Raising up with box Right Hip
subplot(4,1,3)
plot(t2, phase3_LHip_x, t2, phase3_LHip_y,'--', t2, phase3_LHip_z, '.');
legend('x','y','z')
xlabel('Time in seconds');  
ylabel('Degrees');
title('Phase 3');  % Lowering with box Right Hip
subplot(4,1,4)
plot(t5, phase4_LHip_x, t5, phase4_LHip_y,'--', t5, phase4_LHip_z, '.');
legend('x','y','z')
xlabel('Time in seconds'); 
ylabel('Degrees');
title('Phase 4');  % Raising up without box Right Hip

%% ROM calculation

% Right Hip

Rhip_ph1(:,1) = phase1_RHip_x;
Rhip_ph1(:,2) = phase1_RHip_y;
Rhip_ph1(:,3) = phase1_RHip_z;

Rhip_ph2(:,1) = phase2_RHip_x;
Rhip_ph2(:,2) = phase2_RHip_y;
Rhip_ph2(:,3) = phase2_RHip_z;

Rhip_ph3(:,1) = phase3_RHip_x;
Rhip_ph3(:,2) = phase3_RHip_y;
Rhip_ph3(:,3) = phase3_RHip_z;

Rhip_ph4(:,1) = phase4_RHip_x;
Rhip_ph4(:,2) = phase4_RHip_y;
Rhip_ph4(:,3) = phase4_RHip_z;

Rhip_ph1_max = max(Rhip_ph1);
Rhip_ph1_min = min(Rhip_ph1);

Rhip_ph2_max = max(Rhip_ph2);
Rhip_ph2_min = min(Rhip_ph2);

Rhip_ph3_max = max(Rhip_ph3);
Rhip_ph3_min = min(Rhip_ph3);

Rhip_ph4_max = max(Rhip_ph4);
Rhip_ph4_min = min(Rhip_ph4);

for i = 1:length(Rhip_ph1_max)
    Rhip_ph1_ROM(1,i) = Rhip_ph1_max(i)- Rhip_ph1_min(i);
end

for i = 1:length(Rhip_ph1_max)
    Rhip_ph2_ROM(1,i) = Rhip_ph2_max(i)- Rhip_ph2_min(i);
end

for i = 1:length(Rhip_ph1_max)
    Rhip_ph3_ROM(1,i) = Rhip_ph3_max(i)- Rhip_ph3_min(i);
end

for i = 1:length(Rhip_ph1_max)
    Rhip_ph4_ROM(1,i) = Rhip_ph4_max(i)- Rhip_ph4_min(i);
end

% Left Hip

Lhip_ph1(:,1) = phase1_LHip_x;
Lhip_ph1(:,2) = phase1_LHip_y;
Lhip_ph1(:,3) = phase1_LHip_z;

Lhip_ph2(:,1) = phase2_LHip_x;
Lhip_ph2(:,2) = phase2_LHip_y;
Lhip_ph2(:,3) = phase2_LHip_z;

Lhip_ph3(:,1) = phase3_LHip_x;
Lhip_ph3(:,2) = phase3_LHip_y;
Lhip_ph3(:,3) = phase3_LHip_z;

Lhip_ph4(:,1) = phase4_LHip_x;
Lhip_ph4(:,2) = phase4_LHip_y;
Lhip_ph4(:,3) = phase4_LHip_z;

Lhip_ph1_max = max(Lhip_ph1);
Lhip_ph1_min = min(Lhip_ph1);

Lhip_ph2_max = max(Lhip_ph2);
Lhip_ph2_min = min(Lhip_ph2);

Lhip_ph3_max = max(Lhip_ph3);
Lhip_ph3_min = min(Lhip_ph3);

Lhip_ph4_max = max(Lhip_ph4);
Lhip_ph4_min = min(Lhip_ph4);

for i = 1:length(Lhip_ph1_max)
    Lhip_ph1_ROM(1,i) = Lhip_ph1_max(i)- Lhip_ph1_min(i);
end

for i = 1:length(Lhip_ph1_max)
    Lhip_ph2_ROM(1,i) = Lhip_ph2_max(i)- Lhip_ph2_min(i);
end

for i = 1:length(Lhip_ph1_max)
    Lhip_ph3_ROM(1,i) = Lhip_ph3_max(i)- Lhip_ph3_min(i);
end

for i = 1:length(Lhip_ph1_max)
    Lhip_ph4_ROM(1,i) = Lhip_ph4_max(i)- Lhip_ph4_min(i);
end

% Total ROM Hip Angles 

Rhip_total(:,1) = RHipAngles_x;
Rhip_total(:,2) = RHipAngles_y;
Rhip_total(:,3) = RHipAngles_z;

Rhip_total_max = max(Rhip_total);
Rhip_total_min = min(Rhip_total);

for i = 1:length(Rhip_total_max)
    Rhip_total_ROM(1,i) = Rhip_total_max(i)- Rhip_total_min(i);
end

Lhip_total(:,1) = LHipAngles_x;
Lhip_total(:,2) = LHipAngles_y;
Lhip_total(:,3) = LHipAngles_z;

Lhip_total_max = max(Lhip_total);
Lhip_total_min = min(Lhip_total);

for i = 1:length(Lhip_total_max)
    Lhip_total_ROM(1,i) = Lhip_total_max(i)- Lhip_total_min(i);
end

ROM_Right_Hip = table([Rhip_ph1_ROM(1);Rhip_ph2_ROM(1);Rhip_ph3_ROM(1);Rhip_ph4_ROM(1)],[Rhip_ph1_ROM(2);Rhip_ph2_ROM(2);Rhip_ph3_ROM(2);Rhip_ph4_ROM(2)],[Rhip_ph1_ROM(3);Rhip_ph2_ROM(3);Rhip_ph3_ROM(3);Rhip_ph4_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Phase 1','Phase 2','Phase 3', 'Phase 4'})

ROM_Left_Hip = table([Lhip_ph1_ROM(1);Lhip_ph2_ROM(1);Lhip_ph3_ROM(1);Lhip_ph4_ROM(1)],[Lhip_ph1_ROM(2);Lhip_ph2_ROM(2);Lhip_ph3_ROM(2);Lhip_ph4_ROM(2)],[Lhip_ph1_ROM(3);Lhip_ph2_ROM(3);Lhip_ph3_ROM(3);Lhip_ph4_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Phase 1','Phase 2','Phase 3', 'Phase 4'})

TotalROM_Hip = table([Rhip_total_ROM(1); Lhip_total_ROM(1)],[Rhip_total_ROM(2); Lhip_total_ROM(2)],[Rhip_total_ROM(3); Lhip_total_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Right','Left'})



