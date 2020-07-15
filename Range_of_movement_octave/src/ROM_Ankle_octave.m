%% Range of movement of Ankle in Sagittal Lifting

% Code AUTHOR: Yaiza Benito Molpeceres. DATE: January-May 2020.
% Adapted to Octave by Guillermo Asín Prieto

% Segmentation is based on box signal and consists of five phases shown in figure 1:
    % 1: Subject starts at upright position and lowers down to pick the box (idxO1:idxF1)
    % 2: Subject holds and deposits the box (idxF1:idxF2)
    % 3: Subject raises up to upright position and prepares to take the box again (idxF2:idxO3)
    % 4: Subject holds and deposits the box in its original location (idxO3:idxF3)
    % 5: Subject raises up to upright position (idxF3:idxF4)  
    
% ROM is calculated for phases 1, 2, 4 and 5, which are renamed as 1, 2, 3 and 4.clear all % Clear variables

pkg load signal

clear all % Clear variables
close all % Close figures
clc

load('..\tests\data\input\Dinamica44_B.mat')
        
Ts = 1/Fs;                              
t_total = (double(frames)*Ts);
%t_1000 = 0:(t_total/1000):(t_total-(t_total/1000));
t_100 = 1:1:100;

% Normal Ranges of Ankle Motion
for i =1:(length(t_100))
    flex_ankle(1,i) = 45; % flexion
    ext_ankle(1,i) = -20;   % extension
    supination_ankle(1,i) = 20;  % supination
    pronation_ankle(1,i) = -20;  % pronation
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


%% Import model data for Ankle Relative Angles

        RAnkleAngles_x = ModelData.Raw.(ModelOutput{5})(1,:); 
        RAnkleAngles_y = ModelData.Raw.(ModelOutput{5})(2,:);
        RAnkleAngles_z = ModelData.Raw.(ModelOutput{5})(3,:);
        
        LAnkleAngles_x = ModelData.Raw.(ModelOutput{6})(1,:);
        LAnkleAngles_y = ModelData.Raw.(ModelOutput{6})(2,:);
        LAnkleAngles_z = ModelData.Raw.(ModelOutput{6})(3,:);
        
        
%% Ankle signals visualization before going through segmentation

figure(2)
frames = length(position_Z);
% (1:(frames/1000):frames)  si frames>1000
% (1:(frames/1000):frames+(frames/1000)) si frames<1000

plot(t_100,RAnkleAngles_x(1:(frames/100):frames)); 
hold on 
plot(t_100, flex_ankle,'--',t_100, ext_ankle,'--'); 
legend({'x', 'Normal range of flexion-extension'},'Location','best');
%xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
xlabel('% Lifting Cycle');
ylabel('Degrees');
ylim([-30 50]);
title('Right Ankle angles: Flexo-Extension');

figure(3)
plot(t_100, RAnkleAngles_y(1:(frames/100):frames));
hold on 
plot(t_100, supination_ankle,'--',t_100, pronation_ankle,'--'); % ver cual es la positiva
legend({'y', 'Normal range of supination-pronation'},'Location','best');
xlabel('% Lifting Cycle');
ylabel('Degrees');
ylim([-30 30]);
title('Right Ankle angles: Supination-Pronation');

figure(4)
plot(t_100, RAnkleAngles_z(1:(frames/100):frames));
legend({'z'},'Location','best');
xlabel('% Lifting Cycle');
ylabel('Degrees');
%ylim([-70 50]);
title('Right Ankle angles: Internal-External Rotation');

figure(5)
% (1:(frames/100):frames)  si frames>1000
% (1:(frames/100):frames+(frames/100)) si frames<1000

plot(t_100,LAnkleAngles_x(1:(frames/100):frames)); 
hold on 
plot(t_100, flex_ankle,'--',t_100, ext_ankle,'--'); 
legend({'x', 'Normal range of flexion-extension'},'Location','best');
xlabel('% Lifting Cycle');
ylabel('Degrees');
ylim([-30 50]);
title('Left Ankle angles: Flexo-Extension');

figure(6)
plot(t_100, LAnkleAngles_y(1:(frames/100):frames));
hold on 
plot(t_100, supination_ankle,'--',t_100, pronation_ankle,'--'); % ver cual es la positiva
legend({'y', 'Normal range of supination-pronation'},'Location','best');
xlabel('% Lifting Cycle'); 
ylabel('Degrees');
ylim([-30 30]);
title('Left Ankle angles: Supination-Pronation');

figure(7)
plot(t_100, LAnkleAngles_z(1:(frames/100):frames));
legend({'z'},'Location','best');
xlabel('% Lifting Cycle');
ylabel('Degrees');
%ylim([-70 80]);
title('Left Ankle angles: Internal-External Rotation');

%figure(3)
%plot(t ,RAnkleAngles_x, t, RAnkleAngles_y, t, RAnkleAngles_z);
%legend('x','y','z');
%xlabel('time in frames'); % frames totales del trial
%title('Right Ankle angles');

%figure(4)
%plot(t,LAnkleAngles_x, t, LAnkleAngles_y, t, LAnkleAngles_z);   
%legend('x','y','z');
%xlabel('time in frames'); % frames totales del trial
%title('Left Ankle angles');

%% Signal Segmentation

% Right Ankle

phase1_RAnkle_x = RAnkleAngles_x(idxO1:idxF1);  
phase1_RAnkle_y = RAnkleAngles_y(idxO1:idxF1);
phase1_RAnkle_z = RAnkleAngles_z(idxO1:idxF1);

phase2_RAnkle_x = RAnkleAngles_x(idxF1:idxF2); 
phase2_RAnkle_y = RAnkleAngles_y(idxF1:idxF2);
phase2_RAnkle_z = RAnkleAngles_z(idxF1:idxF2);

phase3_RAnkle_x = RAnkleAngles_x(idxO3:idxF3); 
phase3_RAnkle_y = RAnkleAngles_y(idxO3:idxF3);
phase3_RAnkle_z = RAnkleAngles_z(idxO3:idxF3);

phase4_RAnkle_x = RAnkleAngles_x(idxF3:idxF4); 
phase4_RAnkle_y = RAnkleAngles_y(idxF3:idxF4);
phase4_RAnkle_z = RAnkleAngles_z(idxF3:idxF4);

% Left Ankle

phase1_LAnkle_x = LAnkleAngles_x(idxO1:idxF1);
phase1_LAnkle_y = LAnkleAngles_y(idxO1:idxF1);
phase1_LAnkle_z = LAnkleAngles_z(idxO1:idxF1);

phase2_LAnkle_x = LAnkleAngles_x(idxF1:idxF2);
phase2_LAnkle_y = LAnkleAngles_y(idxF1:idxF2);
phase2_LAnkle_z = LAnkleAngles_z(idxF1:idxF2);

phase3_LAnkle_x = LAnkleAngles_x(idxO3:idxF3);
phase3_LAnkle_y = LAnkleAngles_y(idxO3:idxF3);
phase3_LAnkle_z = LAnkleAngles_z(idxO3:idxF3);

phase4_LAnkle_x = LAnkleAngles_x(idxF3:idxF4);
phase4_LAnkle_y = LAnkleAngles_y(idxF3:idxF4);
phase4_LAnkle_z = LAnkleAngles_z(idxF3:idxF4);

%% time vectors definition
t2 = (0:1:length(phase3_RAnkle_x)-1)*Ts;
t3 = (0:1:(length(phase1_RAnkle_x)-1))*Ts;
t4 = (0:1:(length(phase2_RAnkle_x)-1))*Ts;
t5 = (0:1:(length(phase4_RAnkle_x)-1))*Ts;

% Right Ankle figures: MOVEMENT SEGMENTED INTO THE 4 PHASES

figure(8)
subplot(4,1,1)
plot (t3, phase1_RAnkle_x, t3, phase1_RAnkle_y, '--',t3,phase1_RAnkle_z, '.');
legend('x','y','z')
xlabel('Time in seconds');
ylabel('Degrees');
title('RIGHT ANGLES: Phase 1');  % Lowering without box Right Ankle
subplot(4,1,2)
plot(t4, phase2_RAnkle_x,t4,phase2_RAnkle_y, '--',t4,phase2_RAnkle_z, '.');
legend('x','y','z')
xlabel('Time in seconds');
ylabel('Degrees');
title('Phase 2');  %Raising up with box Right Ankle
subplot(4,1,3)
plot(t2, phase3_RAnkle_x, t2, phase3_RAnkle_y,'--', t2, phase3_RAnkle_z, '.');
legend('x','y','z')
xlabel('Time in seconds');
ylabel('Degrees');
title('Phase 3');  % Lowering with box Right Ankle
subplot(4,1,4)
plot(t5, phase4_RAnkle_x, t5, phase4_RAnkle_y,'--', t5, phase4_RAnkle_z, '.');
legend('x','y','z')
xlabel('Time in seconds'); 
ylabel('Degrees');
title('Phase 4');  % Raising up without box Right Ankle

% Left Ankle figures: MOVEMENT SEGMENTED INTO THE 4 PHASES

figure(9)
subplot(4,1,1)
plot (t3, phase1_LAnkle_x, t3, phase1_LAnkle_y, '--',t3,phase1_LAnkle_z, '.');
legend('x','y','z')
xlabel('Time in seconds');
ylabel('Degrees');
title('LEFT ANGLES: Phase 1');  % Lowering without box Right Ankle
subplot(4,1,2)
plot(t4, phase2_LAnkle_x,t4, phase2_LAnkle_y, '--',t4, phase2_LAnkle_z, '.');
legend('x','y','z')
xlabel('Time in seconds');
ylabel('Degrees');
title('Phase 2');  %Raising up with box Right Ankle
subplot(4,1,3)
plot(t2, phase3_LAnkle_x, t2, phase3_LAnkle_y,'--', t2, phase3_LAnkle_z, '.');
legend('x','y','z')
xlabel('Time in seconds');  
ylabel('Degrees');
title('Phase 3');  % Lowering with box Right Ankle
subplot(4,1,4)
plot(t5, phase4_LAnkle_x, t5, phase4_LAnkle_y,'--', t5, phase4_LAnkle_z, '.');
legend('x','y','z')
xlabel('Time in seconds');
ylabel('Degrees');
title('Phase 4');  % Raising up without box Right Ankle

%% ROM calculation

% Right Ankle

RAnkle_ph1(:,1) = phase1_RAnkle_x;
RAnkle_ph1(:,2) = phase1_RAnkle_y;
RAnkle_ph1(:,3) = phase1_RAnkle_z;

RAnkle_ph2(:,1) = phase2_RAnkle_x;
RAnkle_ph2(:,2) = phase2_RAnkle_y;
RAnkle_ph2(:,3) = phase2_RAnkle_z;

RAnkle_ph3(:,1) = phase3_RAnkle_x;
RAnkle_ph3(:,2) = phase3_RAnkle_y;
RAnkle_ph3(:,3) = phase3_RAnkle_z;

RAnkle_ph4(:,1) = phase4_RAnkle_x;
RAnkle_ph4(:,2) = phase4_RAnkle_y;
RAnkle_ph4(:,3) = phase4_RAnkle_z;

RAnkle_ph1_max = max(RAnkle_ph1);
RAnkle_ph1_min = min(RAnkle_ph1);

RAnkle_ph2_max = max(RAnkle_ph2);
RAnkle_ph2_min = min(RAnkle_ph2);

RAnkle_ph3_max = max(RAnkle_ph3);
RAnkle_ph3_min = min(RAnkle_ph3);

RAnkle_ph4_max = max(RAnkle_ph4);
RAnkle_ph4_min = min(RAnkle_ph4);

for i = 1:length(RAnkle_ph1_max)
    RAnkle_ph1_ROM(1,i) = RAnkle_ph1_max(i)- RAnkle_ph1_min(i);
end

for i = 1:length(RAnkle_ph1_max)
    RAnkle_ph2_ROM(1,i) = RAnkle_ph2_max(i)- RAnkle_ph2_min(i);
end

for i = 1:length(RAnkle_ph1_max)
    RAnkle_ph3_ROM(1,i) = RAnkle_ph3_max(i)- RAnkle_ph3_min(i);
end

for i = 1:length(RAnkle_ph1_max)
    RAnkle_ph4_ROM(1,i) = RAnkle_ph4_max(i)- RAnkle_ph4_min(i);
end

% Left Ankle

LAnkle_ph1(:,1) = phase1_LAnkle_x;
LAnkle_ph1(:,2) = phase1_LAnkle_y;
LAnkle_ph1(:,3) = phase1_LAnkle_z;

LAnkle_ph2(:,1) = phase2_LAnkle_x;
LAnkle_ph2(:,2) = phase2_LAnkle_y;
LAnkle_ph2(:,3) = phase2_LAnkle_z;

LAnkle_ph3(:,1) = phase3_LAnkle_x;
LAnkle_ph3(:,2) = phase3_LAnkle_y;
LAnkle_ph3(:,3) = phase3_LAnkle_z;

LAnkle_ph4(:,1) = phase4_LAnkle_x;
LAnkle_ph4(:,2) = phase4_LAnkle_y;
LAnkle_ph4(:,3) = phase4_LAnkle_z;

LAnkle_ph1_max = max(LAnkle_ph1);
LAnkle_ph1_min = min(LAnkle_ph1);

LAnkle_ph2_max = max(LAnkle_ph2);
LAnkle_ph2_min = min(LAnkle_ph2);

LAnkle_ph3_max = max(LAnkle_ph3);
LAnkle_ph3_min = min(LAnkle_ph3);

LAnkle_ph4_max = max(LAnkle_ph4);
LAnkle_ph4_min = min(LAnkle_ph4);

for i = 1:length(LAnkle_ph1_max)
    LAnkle_ph1_ROM(1,i) = LAnkle_ph1_max(i)- LAnkle_ph1_min(i);
end

for i = 1:length(LAnkle_ph1_max)
    LAnkle_ph2_ROM(1,i) = LAnkle_ph2_max(i)- LAnkle_ph2_min(i);
end

for i = 1:length(LAnkle_ph1_max)
    LAnkle_ph3_ROM(1,i) = LAnkle_ph3_max(i)- LAnkle_ph3_min(i);
end

for i = 1:length(LAnkle_ph1_max)
    LAnkle_ph4_ROM(1,i) = LAnkle_ph4_max(i)- LAnkle_ph4_min(i);
end

% Total ROM Ankle Angles

RAnkle_total(:,1) = RAnkleAngles_x;
RAnkle_total(:,2) = RAnkleAngles_y;
RAnkle_total(:,3) = RAnkleAngles_z;

RAnkle_total_max = max(RAnkle_total);
RAnkle_total_min = min(RAnkle_total);

for i = 1:length(RAnkle_total_max)
    RAnkle_total_ROM(1,i) = RAnkle_total_max(i)- RAnkle_total_min(i);
end

LAnkle_total(:,1) = LAnkleAngles_x;
LAnkle_total(:,2) = LAnkleAngles_y;
LAnkle_total(:,3) = LAnkleAngles_z;

LAnkle_total_max = max(LAnkle_total);
LAnkle_total_min = min(LAnkle_total);

for i = 1:length(LAnkle_total_max)
    LAnkle_total_ROM(1,i) = LAnkle_total_max(i)- LAnkle_total_min(i);
end

%ROM_Right_Ankle = table([RAnkle_ph1_ROM(1);RAnkle_ph2_ROM(1);RAnkle_ph3_ROM(1);RAnkle_ph4_ROM(1)],[RAnkle_ph1_ROM(2);RAnkle_ph2_ROM(2);RAnkle_ph3_ROM(2);RAnkle_ph4_ROM(2)],[RAnkle_ph1_ROM(3);RAnkle_ph2_ROM(3);RAnkle_ph3_ROM(3);RAnkle_ph4_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Phase 1','Phase 2','Phase 3', 'Phase 4'})
ROM_Right_Ankle = struct;
ROM_Right_Ankle.Phase1_X=RAnkle_ph1_ROM(1);
ROM_Right_Ankle.Phase1_Y=RAnkle_ph1_ROM(2);
ROM_Right_Ankle.Phase1_Z=RAnkle_ph1_ROM(3);
ROM_Right_Ankle.Phase2_X=RAnkle_ph2_ROM(1);
ROM_Right_Ankle.Phase2_Y=RAnkle_ph2_ROM(2);
ROM_Right_Ankle.Phase2_Z=RAnkle_ph2_ROM(3);
ROM_Right_Ankle.Phase3_X=RAnkle_ph3_ROM(1);
ROM_Right_Ankle.Phase3_Y=RAnkle_ph3_ROM(2);
ROM_Right_Ankle.Phase3_Z=RAnkle_ph3_ROM(3);
ROM_Right_Ankle.Phase4_X=RAnkle_ph4_ROM(1);
ROM_Right_Ankle.Phase4_Y=RAnkle_ph4_ROM(2);
ROM_Right_Ankle.Phase4_Z=RAnkle_ph4_ROM(3);   

ROM_Right_Ankle

%ROM_Left_Ankle = table([LAnkle_ph1_ROM(1);LAnkle_ph2_ROM(1);LAnkle_ph3_ROM(1);LAnkle_ph4_ROM(1)],[LAnkle_ph1_ROM(2);LAnkle_ph2_ROM(2);LAnkle_ph3_ROM(2);LAnkle_ph4_ROM(2)],[LAnkle_ph1_ROM(3);LAnkle_ph2_ROM(3);LAnkle_ph3_ROM(3);LAnkle_ph4_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Phase 1','Phase 2','Phase 3', 'Phase 4'})
ROM_Left_Ankle = struct;
ROM_Left_Ankle.Phase1_X=LAnkle_ph1_ROM(1);
ROM_Left_Ankle.Phase1_Y=LAnkle_ph1_ROM(2);
ROM_Left_Ankle.Phase1_Z=LAnkle_ph1_ROM(3);
ROM_Left_Ankle.Phase2_X=LAnkle_ph2_ROM(1);
ROM_Left_Ankle.Phase2_Y=LAnkle_ph2_ROM(2);
ROM_Left_Ankle.Phase2_Z=LAnkle_ph2_ROM(3);
ROM_Left_Ankle.Phase3_X=LAnkle_ph3_ROM(1);
ROM_Left_Ankle.Phase3_Y=LAnkle_ph3_ROM(2);
ROM_Left_Ankle.Phase3_Z=LAnkle_ph3_ROM(3);
ROM_Left_Ankle.Phase4_X=LAnkle_ph4_ROM(1);
ROM_Left_Ankle.Phase4_Y=LAnkle_ph4_ROM(2);
ROM_Left_Ankle.Phase4_Z=LAnkle_ph4_ROM(3);

ROM_Left_Ankle

%TotalROM_Ankle = table([RAnkle_total_ROM(1); LAnkle_total_ROM(1)],[RAnkle_total_ROM(2); LAnkle_total_ROM(2)],[RAnkle_total_ROM(3); LAnkle_total_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Right','Left'})
TotalROM_Ankle = struct;
TotalROM_Ankle.Right_X=RAnkle_total_ROM(1);
TotalROM_Ankle.Right_Y=RAnkle_total_ROM(2);
TotalROM_Ankle.Right_Z=RAnkle_total_ROM(3);
TotalROM_Ankle.Left_X=LAnkle_total_ROM(1);
TotalROM_Ankle.Left_Y=LAnkle_total_ROM(2);
TotalROM_Ankle.Left_Z=LAnkle_total_ROM(3);

TotalROM_Ankle
