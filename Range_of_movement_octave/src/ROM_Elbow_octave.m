%% Range of movement of Elbow in Sagittal Lifting

% Code AUTHOR: Yaiza Benito Molpeceres. DATE: January-May 2020.
% Adapted to Octave by Guillermo Asín Prieto

% Segmentation is based on box signal and consists of five phases shown in figure 1:
    % 1: Subject starts at upright position and lowers down to pick the box (idxO1:idxF1)
    % 2: Subject holds and deposits the box (idxF1:idxF2)
    % 3: Subject raises up to upright position and prepares to take the box again (idxF2:idxO3)
    % 4: Subject holds and deposits the box in its original location (idxO3:idxF3)
    % 5: Subject raises up to upright position (idxF3:idxF4)  
    
% ROM is calculated for phases 1, 2, 4 and 5, which are renamed as 1, 2, 3 and 4.

pkg load signal

clear all % Clear variables
close all % Close figures
clc

load('..\tests\data\input\dinamica44_B.mat')
            
Ts = 1/Fs;
t_total = (double(frames)*Ts);
%t_1000 = 0:(t_total/1000):(t_total-(t_total/1000));
t_100 = 1:1:100;

% Normal Ranges of Elbow Motion
for i =1:(length(t_100))
    flex_elbow(1,i) = 130; % flexion o 180 V
    ext_elbow(1,i) = -45;   % extension o 60 V
    supination_elbow(1,i) = 20;  % supination
    pronation_elbow(1,i) = -20;  % pronation
end

%% Events of interest in vertical box trajectory identification
 
position_Z = SIDEBOX3(:,3)';
 
%[pks, locs] = findpeaks(position_Z,'minPeakProminence',20);
[pks, locs] = findpeaks(position_Z,'DoubleSided');   
%locs=locs(pks>80);                 
%pks=pks(pks>80);
for count=1:1:length(pks)
    if count==1
        signal_aux_L=position_Z(1:locs(count));
        signal_aux_R=position_Z(locs(count):locs(count+1));
    else if count==length(pks)
            signal_aux_L=position_Z(locs(count-1):locs(count));
            signal_aux_R=position_Z(locs(count):end);
        else
            signal_aux_L=position_Z(locs(count-1):locs(count));
            signal_aux_R=position_Z(locs(count):locs(count+1));
        end
    end
    try
      min_L_aux=findpeaks(1.01*max(signal_aux_L)-signal_aux_L,'DoubleSided');
      min_L=min(1.01*max(signal_aux_L)-min_L_aux);
      if isempty(min_L)
        min_L=min(signal_aux_L);
      end
    catch
      min_L=min(signal_aux_L);
    end
    try
      min_R_aux=findpeaks(1.01*max(signal_aux_R)-signal_aux_R,'DoubleSided');
      min_R=min(1.01*max(signal_aux_R)-min_R_aux);
      if isempty(min_R)
        min_R=min(signal_aux_R);
      end
    catch
      min_R=min(signal_aux_R);
    end
    real_min=min(min_L,min_R);
    prominence(count)=pks(count)-real_min;
end
pks=pks(prominence>20);    
locs=locs(prominence>20);


t = (0:1:(length(position_Z)-1))*Ts;

%TF1 = islocalmin(position_Z, 'FlatSelection','first');
[TF1,locs_TF1] = findpeaks(1.01*max(position_Z(1:locs(1)))-position_Z(1:locs(1)),'DoubleSided');
TF1=1.01*max(position_Z(1:locs(1)))-TF1;

TF1_aux=zeros(1,length(TF1));
for count=1:1:length(TF1)-1
  if TF1(count+1)==TF1(count)
    TF1_aux(count+1)=1;
  end
end
TF1(find(TF1_aux))=[];
locs_TF1(find(TF1_aux))=[]; 


for count_aux=1:1:length(pks)-1
  %TF1 = islocalmin(position_Z, 'FlatSelection','first');
  [TF1_2,locs_TF1_2] = findpeaks(1.01*max(position_Z(locs(count_aux):locs(count_aux+1)))-position_Z(locs(count_aux):locs(count_aux+1)),'DoubleSided');
  TF1_2=1.01*max(position_Z(locs(count_aux):locs(count_aux+1)))-TF1_2;

  TF1_aux_2=zeros(1,length(TF1_2));
  for count=1:1:length(TF1_2)-1
    if TF1_2(count+1)==TF1_2(count)
      TF1_aux_2(count+1)=1;
    end
  end
  TF1_2(find(TF1_aux_2))=[];
  locs_TF1_2(find(TF1_aux_2))=[];
  TF1 = [TF1,TF1_2];
  locs_TF1 = [locs_TF1,locs_TF1_2+locs(count_aux)];
endfor


%TF1 = islocalmin(position_Z, 'FlatSelection','first');
[TF1_2,locs_TF1_2] = findpeaks(1.01*max(position_Z(locs(end):end))-position_Z(locs(end):end),'DoubleSided');
TF1_2=1.01*max(position_Z(locs(end):end))-TF1_2;

TF1_aux_2=zeros(1,length(TF1_2));
for count=1:1:length(TF1_2)-1
  if TF1_2(count+1)==TF1_2(count)
    TF1_aux_2(count+1)=1;
  end
end
TF1_2(find(TF1_aux_2))=[];
locs_TF1_2(find(TF1_aux_2))=[]; 

TF1 = [TF1,TF1_2];
locs_TF1 = [locs_TF1,locs_TF1_2+locs(end)];

for count=1:1:length(TF1)-1
    for count_aux=1:1:length(pks)
        if TF1(count)==pks(count_aux)
          TF1_aux(count)=1;
        end
    end
end

TF1(find(TF1_aux))=[];
locs_TF1(find(TF1_aux))=[];


%for i = 1:frames                               % solo para dinamica 03
 %   TF2(1,i) = position_Z(1,i) > 0;
%end
%idx = find(TF1);
idx = locs_TF1;

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
%plot(t,position_Z,'r*','MarkerIndices',idxO1);
plot(t(idxO1),position_Z(idxO1),'r*');

%plot(t,position_Z,'r*','MarkerIndices',idxF1);
plot(t(idxF1),position_Z(idxF1),'r*');

%plot(t,position_Z,'r*','MarkerIndices',idxF2);
plot(t(idxF2),position_Z(idxF2),'r*');

%plot(t,position_Z,'r*','MarkerIndices',idxO3);
plot(t(idxO3),position_Z(idxO3),'r*');

%plot(t,position_Z,'r*','MarkerIndices',idxF3);
plot(t(idxF3),position_Z(idxF3),'r*');

%plot(t,position_Z,'r*','MarkerIndices',idxF4);
plot(t(idxF4),position_Z(idxF4),'r*');

hold off
xlabel('Time in seconds');
ylabel('Millimetres');
title('Box signal segmentation in 4 phases'); % or title: Sagittal Lifting Segmentation


%% Import model data for Elbow Relative Angles

        RElbowAngles_x = ModelData.Raw.(ModelOutput{7})(1,:); 
        RElbowAngles_y = ModelData.Raw.(ModelOutput{7})(2,:);
        RElbowAngles_z = ModelData.Raw.(ModelOutput{7})(3,:);
        
        LElbowAngles_x = ModelData.Raw.(ModelOutput{8})(1,:);
        LElbowAngles_y = ModelData.Raw.(ModelOutput{8})(2,:);
        LElbowAngles_z = ModelData.Raw.(ModelOutput{8})(3,:);
        

%% Elbow signals visualization before going through segmentation

figure(2)
frames = length(position_Z);
% (1:(frames/100):frames)  si frames>1000
% (1:(frames/100):frames+(frames/100)) si frames<1000

plot(t_100,RElbowAngles_x(1:(frames/100):frames)); 
hold on 
plot(t_100, flex_elbow,'--',t_100, ext_elbow,'--'); 
legend({'x', 'Normal range of flexion-extension'},'Location','best');
xlabel('% Lifting Cycle'); 
ylabel('Degrees');
%ylim([-100 160]);
title('Right Elbow angles: Flexo-Extension');

figure(3)
plot(t_100, RElbowAngles_y(1:(frames/100):frames));
legend({'y'},'Location','best');
xlabel('% Lifting Cycle');  
ylabel('Degrees');
%ylim([-70 50]);
title('Right Elbow angles: Supination-Pronation');

figure(4)
plot(t_100, RElbowAngles_z(1:(frames/100):frames));
legend({'z'},'Location','best');
xlabel('% Lifting Cycle'); 
ylabel('Degrees');
%ylim([-70 50]);
title('Right Elbow angles: Internal-External Rotation');

figure(5)
% (1:(frames/100):frames)  si frames>1000
% (1:(frames/100):frames+(frames/100)) si frames<1000

plot(t_100,LElbowAngles_x(1:(frames/100):frames)); 
hold on 
plot(t_100, flex_elbow,'--',t_100, ext_elbow,'--'); 
legend({'x', 'Normal range of flexion-extension'},'Location','best');
xlabel('% Lifting Cycle'); 
ylabel('Degrees');
%ylim([-100 160]);
title('Left Elbow angles: Flexo-Extension');

figure(6)
plot(t_100, LElbowAngles_y(1:(frames/100):frames));
legend({'y'},'Location','best');
xlabel('% Lifting Cycle'); 
ylabel('Degrees');
%ylim([-70 50]);
title('Left Elbow angles: Supination-Pronation');

figure(7)
plot(t_100, LElbowAngles_z(1:(frames/100):frames));
legend({'z'},'Location','best');
xlabel('% Lifting Cycle'); 
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

phase3_RElbow_x = RElbowAngles_x(idxO3:idxF3); 
phase3_RElbow_y = RElbowAngles_y(idxO3:idxF3);
phase3_RElbow_z = RElbowAngles_z(idxO3:idxF3);

phase4_RElbow_x = RElbowAngles_x(idxF3:idxF4); 
phase4_RElbow_y = RElbowAngles_y(idxF3:idxF4);
phase4_RElbow_z = RElbowAngles_z(idxF3:idxF4);

% Left Elbow 

phase1_LElbow_x = LElbowAngles_x(idxO1:idxF1);
phase1_LElbow_y = LElbowAngles_y(idxO1:idxF1);
phase1_LElbow_z = LElbowAngles_z(idxO1:idxF1);

phase2_LElbow_x = LElbowAngles_x(idxF1:idxF2);
phase2_LElbow_y = LElbowAngles_y(idxF1:idxF2);
phase2_LElbow_z = LElbowAngles_z(idxF1:idxF2);

phase3_LElbow_x = LElbowAngles_x(idxO3:idxF3);
phase3_LElbow_y = LElbowAngles_y(idxO3:idxF3);
phase3_LElbow_z = LElbowAngles_z(idxO3:idxF3);

phase4_LElbow_x = LElbowAngles_x(idxF3:idxF4);
phase4_LElbow_y = LElbowAngles_y(idxF3:idxF4);
phase4_LElbow_z = LElbowAngles_z(idxF3:idxF4);

%% time vectors definition

t2 = (0:1:(length(phase3_RElbow_x)-1))*Ts;
t3 = (0:1:(length(phase1_RElbow_x)-1))*Ts;
t4 = (0:1:(length(phase2_RElbow_x)-1))*Ts;
t5 = (0:1:(length(phase4_RElbow_x)-1))*Ts;

% Right Elbow figures: MOVEMENT SEGMENTED INTO THE 4 PHASES

figure(8)
subplot(4,1,1)
plot (t3, phase1_RElbow_x, t3, phase1_RElbow_y, '--',t3,phase1_RElbow_z, '.');
legend('x','y','z')
xlabel('Time in seconds'); 
ylabel('Degrees');
title('Phase 1');  % Lowering without box Right Elbow
subplot(4,1,2)
plot(t4, phase2_RElbow_x,t4,phase2_RElbow_y, '--',t4,phase2_RElbow_z, '.');
legend('x','y','z')
xlabel('Time in seconds'); 
ylabel('Degrees');
title('Phase 2');  %Raising up with box Right Elbow
subplot(4,1,3)
plot(t2, phase3_RElbow_x, t2, phase3_RElbow_y,'--', t2, phase3_RElbow_z, '.');
legend('x','y','z')
xlabel('Time in seconds'); 
ylabel('Degrees');
title('Phase 3');  % Lowering with box Right Elbow
subplot(4,1,4)
plot(t5, phase4_RElbow_x, t5, phase4_RElbow_y,'--', t5, phase4_RElbow_z, '.');
legend('x','y','z')
xlabel('Time in seconds'); 
ylabel('Degrees');
title('Phase 4');  % Raising up without box Right Elbow

% Left Elbow figures: MOVEMENT SEGMENTED INTO THE 4 PHASES

figure(9)
subplot(4,1,1)
plot (t3, phase1_LElbow_x, t3, phase1_LElbow_y, '--',t3,phase1_LElbow_z, '.');
legend('x','y','z')
xlabel('Time in seconds'); 
ylabel('Degrees');
title('Phase 1');  % Lowering without box Right Elbow
subplot(4,1,2)
plot(t4, phase2_LElbow_x,t4, phase2_LElbow_y, '--',t4, phase2_LElbow_z, '.');
legend('x','y','z')
xlabel('Time in seconds');  
ylabel('Degrees');
title('Phase 2');  %Raising up with box Right Elbow
subplot(4,1,3)
plot(t2, phase3_LElbow_x, t2, phase3_LElbow_y,'--', t2, phase3_LElbow_z, '.');
legend('x','y','z')
xlabel('Time in seconds'); 
ylabel('Degrees');
title('Phase 3');  % Lowering with box Right Elbow
subplot(4,1,4)
plot(t5, phase4_LElbow_x, t5, phase4_LElbow_y,'--', t5, phase4_LElbow_z, '.');
legend('x','y','z')
xlabel('Time in seconds'); 
ylabel('Degrees');
title('Phase 4');  % Raising up without box Right Elbow

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

RElbow_ph1_max = max(RElbow_ph1);
RElbow_ph1_min = min(RElbow_ph1);

RElbow_ph2_max = max(RElbow_ph2);
RElbow_ph2_min = min(RElbow_ph2);

RElbow_ph3_max = max(RElbow_ph3);
RElbow_ph3_min = min(RElbow_ph3);

RElbow_ph4_max = max(RElbow_ph4);
RElbow_ph4_min = min(RElbow_ph4);

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

LElbow_ph1_max = max(LElbow_ph1);
LElbow_ph1_min = min(LElbow_ph1);

LElbow_ph2_max = max(LElbow_ph2);
LElbow_ph2_min = min(LElbow_ph2);

LElbow_ph3_max = max(LElbow_ph3);
LElbow_ph3_min = min(LElbow_ph3);

LElbow_ph4_max = max(LElbow_ph4);
LElbow_ph4_min = min(LElbow_ph4);

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

%ROM_Right_Elbow = table([RElbow_ph1_ROM(1);RElbow_ph2_ROM(1);RElbow_ph3_ROM(1);RElbow_ph4_ROM(1)],[RElbow_ph1_ROM(2);RElbow_ph2_ROM(2);RElbow_ph3_ROM(2);RElbow_ph4_ROM(2)],[RElbow_ph1_ROM(3);RElbow_ph2_ROM(3);RElbow_ph3_ROM(3);RElbow_ph4_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Phase 1','Phase 2','Phase 3', 'Phase 4'})
ROM_Right_Elbow = struct;
ROM_Right_Elbow.Phase1_X=RElbow_ph1_ROM(1);
ROM_Right_Elbow.Phase1_Y=RElbow_ph1_ROM(2);
ROM_Right_Elbow.Phase1_Z=RElbow_ph1_ROM(3);
ROM_Right_Elbow.Phase2_X=RElbow_ph2_ROM(1);
ROM_Right_Elbow.Phase2_Y=RElbow_ph2_ROM(2);
ROM_Right_Elbow.Phase2_Z=RElbow_ph2_ROM(3);
ROM_Right_Elbow.Phase3_X=RElbow_ph3_ROM(1);
ROM_Right_Elbow.Phase3_Y=RElbow_ph3_ROM(2);
ROM_Right_Elbow.Phase3_Z=RElbow_ph3_ROM(3);
ROM_Right_Elbow.Phase4_X=RElbow_ph4_ROM(1);
ROM_Right_Elbow.Phase4_Y=RElbow_ph4_ROM(2);
ROM_Right_Elbow.Phase4_Z=RElbow_ph4_ROM(3);   

ROM_Right_Elbow

%ROM_Left_Elbow = table([LElbow_ph1_ROM(1);LElbow_ph2_ROM(1);LElbow_ph3_ROM(1);LElbow_ph4_ROM(1)],[LElbow_ph1_ROM(2);LElbow_ph2_ROM(2);LElbow_ph3_ROM(2);LElbow_ph4_ROM(2)],[LElbow_ph1_ROM(3);LElbow_ph2_ROM(3);LElbow_ph3_ROM(3);LElbow_ph4_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Phase 1','Phase 2','Phase 3', 'Phase 4'})
ROM_Left_Elbow = struct;
ROM_Left_Elbow.Phase1_X=LElbow_ph1_ROM(1);
ROM_Left_Elbow.Phase1_Y=LElbow_ph1_ROM(2);
ROM_Left_Elbow.Phase1_Z=LElbow_ph1_ROM(3);
ROM_Left_Elbow.Phase2_X=LElbow_ph2_ROM(1);
ROM_Left_Elbow.Phase2_Y=LElbow_ph2_ROM(2);
ROM_Left_Elbow.Phase2_Z=LElbow_ph2_ROM(3);
ROM_Left_Elbow.Phase3_X=LElbow_ph3_ROM(1);
ROM_Left_Elbow.Phase3_Y=LElbow_ph3_ROM(2);
ROM_Left_Elbow.Phase3_Z=LElbow_ph3_ROM(3);
ROM_Left_Elbow.Phase4_X=LElbow_ph4_ROM(1);
ROM_Left_Elbow.Phase4_Y=LElbow_ph4_ROM(2);
ROM_Left_Elbow.Phase4_Z=LElbow_ph4_ROM(3);

ROM_Left_Elbow

%TotalROM_Elbow = table([RElbow_total_ROM(1); LElbow_total_ROM(1)],[RElbow_total_ROM(2); LElbow_total_ROM(2)],[RElbow_total_ROM(3); LElbow_total_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Right','Left'})
TotalROM_Elbow = struct;
TotalROM_Elbow.Right_X=RElbow_total_ROM(1);
TotalROM_Elbow.Right_Y=RElbow_total_ROM(2);
TotalROM_Elbow.Right_Z=RElbow_total_ROM(3);
TotalROM_Elbow.Left_X=LElbow_total_ROM(1);
TotalROM_Elbow.Left_Y=LElbow_total_ROM(2);
TotalROM_Elbow.Left_Z=LElbow_total_ROM(3);

TotalROM_Elbow