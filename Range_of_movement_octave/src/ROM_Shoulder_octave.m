%% Range of movement of Shoulder in Sagittal Lifting

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

load('..\tests\data\input\Dinamica44_B.mat')

Ts = 1/Fs;
t_total = (double(frames)*Ts);
%t_1000 = 0:(t_total/1000):(t_total-(t_total/1000));
t_100 = 1:1:100;

% Normal Ranges of Shoulder Motion
for i =1:(length(t_100))
    flex_shoulder(1,i) = 130; % flexion
    ext_shoulder(1,i) = -30;   % extension
    abd_shoulder(1,i) = -50;   % abduction
    add_shoulder(1,i) = 30;   % adduction
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


%% Import model data for Shoulder Relative Angles

        RShoulderAngles_x = ModelData.Raw.(ModelOutput{9})(1,:); 
        RShoulderAngles_y = ModelData.Raw.(ModelOutput{9})(2,:);
        RShoulderAngles_z = ModelData.Raw.(ModelOutput{9})(3,:);
        
        LShoulderAngles_x = ModelData.Raw.(ModelOutput{10})(1,:);
        LShoulderAngles_y = ModelData.Raw.(ModelOutput{10})(2,:);
        LShoulderAngles_z = ModelData.Raw.(ModelOutput{10})(3,:);
        
        
%% Shoulder signals visualization before going through segmentation

figure(2)
frames = length(position_Z);
% (1:(frames/100):frames)  si frames>1000
% (1:(frames/100):frames+(frames/100)) si frames<1000

plot(t_100,RShoulderAngles_x(1:(frames/100):frames)); 
hold on 
plot(t_100, flex_shoulder,'--',t_100, ext_shoulder,'--'); 
legend({'x', 'Normal range of flexion-extension'},'Location','best');
xlabel('% Lifting Cycle');
ylabel('Degrees');
%ylim([-100 160]);
title('Right Shoulder angles: Flexo-Extension');

figure(3)
plot(t_100, RShoulderAngles_y(1:(frames/100):frames));
hold on 
plot(t_100, abd_shoulder,'--',t_100, add_shoulder,'--'); 
legend({'y', 'Normal range of abduction-adduction'},'Location','best');
xlabel('% Lifting Cycle');
ylabel('Degrees');
ylim([-60 100]);
title('Right Shoulder angles: Abduction-Adduction');

figure(4)
plot(t_100, RShoulderAngles_z(1:(frames/100):frames));
legend({'z'},'Location','best');
xlabel('% Lifting Cycle');
ylabel('Degrees');
%ylim([-70 50]);
title('Right Shoulder angles: Internal-External Rotation');

figure(5)
% (1:(frames/100):frames)  si frames>1000
% (1:(frames/100):frames+(frames/100)) si frames<1000

plot(t_100,LShoulderAngles_x(1:(frames/100):frames)); 
hold on 
plot(t_100, flex_shoulder,'--',t_100, ext_shoulder,'--'); 
legend({'x', 'Normal range of flexion-extension'},'Location','best');
xlabel('% Lifting Cycle');
ylabel('Degrees');
%ylim([-100 160]);
title('Left Shoulder angles: Flexo-Extension');

figure(6)
plot(t_100, LShoulderAngles_y(1:(frames/100):frames));
hold on 
plot(t_100, abd_shoulder,'--',t_100, add_shoulder,'--'); 
legend({'y', 'Normal range of abduction-adduction'},'Location','best');
xlabel('% Lifting Cycle');
ylabel('Degrees');
ylim([-60 100]);
title('Left Shoulder angles: Abduction-Adduction');

figure(7)
plot(t_100, LShoulderAngles_z(1:(frames/100):frames));
legend({'z', 'Normal range of internal-external rotation'},'Location','best');
xlabel('% Lifting Cycle'); 
ylabel('Degrees');
%ylim([-70 80]);
title('Left Shoulder angles: Internal-External Rotation');

%figure(3)
%plot(t ,RShoulderAngles_x, t, RShoulderAngles_y, t, RShoulderAngles_z);
%legend('x','y','z');
%xlabel('time in frames'); % frames totales del trial
%title('Right Shoulder angles');

%figure(4)
%plot(t,LShoulderAngles_x, t, LShoulderAngles_y, t, LShoulderAngles_z);   
%legend('x','y','z');
%xlabel('time in frames'); % frames totales del trial
%title('Left Shoulder angles');

%% Signal Segmentation

% Right Shoulder

phase1_RShoulder_x = RShoulderAngles_x(idxO1:idxF1);  
phase1_RShoulder_y = RShoulderAngles_y(idxO1:idxF1);
phase1_RShoulder_z = RShoulderAngles_z(idxO1:idxF1);

phase2_RShoulder_x = RShoulderAngles_x(idxF1:idxF2); 
phase2_RShoulder_y = RShoulderAngles_y(idxF1:idxF2);
phase2_RShoulder_z = RShoulderAngles_z(idxF1:idxF2);

phase3_RShoulder_x = RShoulderAngles_x(idxO3:idxF3); 
phase3_RShoulder_y = RShoulderAngles_y(idxO3:idxF3);
phase3_RShoulder_z = RShoulderAngles_z(idxO3:idxF3);

phase4_RShoulder_x = RShoulderAngles_x(idxF3:idxF4); 
phase4_RShoulder_y = RShoulderAngles_y(idxF3:idxF4);
phase4_RShoulder_z = RShoulderAngles_z(idxF3:idxF4);

% Left Shoulder

phase1_LShoulder_x = LShoulderAngles_x(idxO1:idxF1);
phase1_LShoulder_y = LShoulderAngles_y(idxO1:idxF1);
phase1_LShoulder_z = LShoulderAngles_z(idxO1:idxF1);

phase2_LShoulder_x = LShoulderAngles_x(idxF1:idxF2);
phase2_LShoulder_y = LShoulderAngles_y(idxF1:idxF2);
phase2_LShoulder_z = LShoulderAngles_z(idxF1:idxF2);

phase3_LShoulder_x = LShoulderAngles_x(idxO3:idxF3);
phase3_LShoulder_y = LShoulderAngles_y(idxO3:idxF3);
phase3_LShoulder_z = LShoulderAngles_z(idxO3:idxF3);

phase4_LShoulder_x = LShoulderAngles_x(idxF3:idxF4);
phase4_LShoulder_y = LShoulderAngles_y(idxF3:idxF4);
phase4_LShoulder_z = LShoulderAngles_z(idxF3:idxF4);

%% time vectors definition

t2 = (0:1:(length(phase3_RShoulder_x)-1))*Ts;
t3 = (0:1:(length(phase1_RShoulder_x)-1))*Ts;
t4 = (0:1:(length(phase2_RShoulder_x)-1))*Ts;
t5 = (0:1:(length(phase4_RShoulder_x)-1))*Ts;

% Right Shoulder figures: MOVEMENT SEGMENTED INTO THE 4 PHASES

figure(8)
subplot(4,1,1)
plot (t3, phase1_RShoulder_x, t3, phase1_RShoulder_y, '--',t3,phase1_RShoulder_z, '.');
legend('x','y','z')
xlabel('Time in seconds'); 
ylabel('Degrees');
title('RIGHT ANGLES: Phase 1');  % Lowering without box Right Shoulder
subplot(4,1,2)
plot(t4, phase2_RShoulder_x,t4,phase2_RShoulder_y, '--',t4,phase2_RShoulder_z, '.');
legend('x','y','z')
xlabel('Time in seconds'); 
ylabel('Degrees');
title('Phase 2');  %Raising up with box Right Shoulder
subplot(4,1,3)
plot(t2, phase3_RShoulder_x, t2, phase3_RShoulder_y,'--', t2, phase3_RShoulder_z, '.');
legend('x','y','z')
xlabel('Time in seconds');  
ylabel('Degrees');
title('Phase 3');  % Lowering with box Right Shoulder
subplot(4,1,4)
plot(t5, phase4_RShoulder_x, t5, phase4_RShoulder_y,'--', t5, phase4_RShoulder_z, '.');
legend('x','y','z')
xlabel('Time in seconds'); 
ylabel('Degrees');
title('Phase 4');  % Raising up without box Right Shoulder

% Left Shoulder figures: MOVEMENT SEGMENTED INTO THE 4 PHASES

figure(9)
subplot(4,1,1)
plot (t3, phase1_LShoulder_x, t3, phase1_LShoulder_y, '--',t3,phase1_LShoulder_z, '.');
legend('x','y','z')
xlabel('Time in seconds'); 
ylabel('Degrees');
title('LEFT ANGLES: Phase 1');  % Lowering without box Right Shoulder
subplot(4,1,2)
plot(t4, phase2_LShoulder_x,t4, phase2_LShoulder_y, '--',t4, phase2_LShoulder_z, '.');
legend('x','y','z')
xlabel('Time in seconds'); 
ylabel('Degrees');
title('Phase 2');  %Raising up with box Right Shoulder
subplot(4,1,3)
plot(t2, phase3_LShoulder_x, t2, phase3_LShoulder_y,'--', t2, phase3_LShoulder_z, '.');
legend('x','y','z')
xlabel('Time in seconds'); 
ylabel('Degrees');
title('Phase 3');  % Lowering with box Right Shoulder
subplot(4,1,4)
plot(t5, phase4_LShoulder_x, t5, phase4_LShoulder_y,'--', t5, phase4_LShoulder_z, '.');
legend('x','y','z')
xlabel('Time in seconds');  
ylabel('Degrees');
title('Phase 4');  % Raising up without box Right Shoulder

%% ROM calculation

% Right Shoulder

RShoulder_ph1(:,1) = phase1_RShoulder_x;
RShoulder_ph1(:,2) = phase1_RShoulder_y;
RShoulder_ph1(:,3) = phase1_RShoulder_z;

RShoulder_ph2(:,1) = phase2_RShoulder_x;
RShoulder_ph2(:,2) = phase2_RShoulder_y;
RShoulder_ph2(:,3) = phase2_RShoulder_z;

RShoulder_ph3(:,1) = phase3_RShoulder_x;
RShoulder_ph3(:,2) = phase3_RShoulder_y;
RShoulder_ph3(:,3) = phase3_RShoulder_z;

RShoulder_ph4(:,1) = phase4_RShoulder_x;
RShoulder_ph4(:,2) = phase4_RShoulder_y;
RShoulder_ph4(:,3) = phase4_RShoulder_z;

RShoulder_ph1_max = max(RShoulder_ph1);
RShoulder_ph1_min = min(RShoulder_ph1);

RShoulder_ph2_max = max(RShoulder_ph2);
RShoulder_ph2_min = min(RShoulder_ph2);

RShoulder_ph3_max = max(RShoulder_ph3);
RShoulder_ph3_min = min(RShoulder_ph3);

RShoulder_ph4_max = max(RShoulder_ph4);
RShoulder_ph4_min = min(RShoulder_ph4);

for i = 1:length(RShoulder_ph1_max)
    RShoulder_ph1_ROM(1,i) = RShoulder_ph1_max(i)- RShoulder_ph1_min(i);
end

for i = 1:length(RShoulder_ph1_max)
    RShoulder_ph2_ROM(1,i) = RShoulder_ph2_max(i)- RShoulder_ph2_min(i);
end

for i = 1:length(RShoulder_ph1_max)
    RShoulder_ph3_ROM(1,i) = RShoulder_ph3_max(i)- RShoulder_ph3_min(i);
end

for i = 1:length(RShoulder_ph1_max)
    RShoulder_ph4_ROM(1,i) = RShoulder_ph4_max(i)- RShoulder_ph4_min(i);
end

% Left Shoulder

LShoulder_ph1(:,1) = phase1_LShoulder_x;
LShoulder_ph1(:,2) = phase1_LShoulder_y;
LShoulder_ph1(:,3) = phase1_LShoulder_z;

LShoulder_ph2(:,1) = phase2_LShoulder_x;
LShoulder_ph2(:,2) = phase2_LShoulder_y;
LShoulder_ph2(:,3) = phase2_LShoulder_z;

LShoulder_ph3(:,1) = phase3_LShoulder_x;
LShoulder_ph3(:,2) = phase3_LShoulder_y;
LShoulder_ph3(:,3) = phase3_LShoulder_z;

LShoulder_ph4(:,1) = phase4_LShoulder_x;
LShoulder_ph4(:,2) = phase4_LShoulder_y;
LShoulder_ph4(:,3) = phase4_LShoulder_z;

LShoulder_ph1_max = max(LShoulder_ph1);
LShoulder_ph1_min = min(LShoulder_ph1);

LShoulder_ph2_max = max(LShoulder_ph2);
LShoulder_ph2_min = min(LShoulder_ph2);

LShoulder_ph3_max = max(LShoulder_ph3);
LShoulder_ph3_min = min(LShoulder_ph3);

LShoulder_ph4_max = max(LShoulder_ph4);
LShoulder_ph4_min = min(LShoulder_ph4);

for i = 1:length(RShoulder_ph1_max)
    LShoulder_ph1_ROM(1,i) = LShoulder_ph1_max(i)- LShoulder_ph1_min(i);
end

for i = 1:length(RShoulder_ph1_max)
    LShoulder_ph2_ROM(1,i) = LShoulder_ph2_max(i)- LShoulder_ph2_min(i);
end

for i = 1:length(RShoulder_ph1_max)
    LShoulder_ph3_ROM(1,i) = LShoulder_ph3_max(i)- LShoulder_ph3_min(i);
end

for i = 1:length(RShoulder_ph1_max)
    LShoulder_ph4_ROM(1,i) = LShoulder_ph4_max(i)- LShoulder_ph4_min(i);
end

% Total ROM Shoulder Angles

RShoulder_total(:,1) = RShoulderAngles_x;
RShoulder_total(:,2) = RShoulderAngles_y;
RShoulder_total(:,3) = RShoulderAngles_z;

RShoulder_total_max = max(RShoulder_total);
RShoulder_total_min = min(RShoulder_total);

for i = 1:length(RShoulder_total_max)
    RShoulder_total_ROM(1,i) = RShoulder_total_max(i)- RShoulder_total_min(i);
end

LShoulder_total(:,1) = LShoulderAngles_x;
LShoulder_total(:,2) = LShoulderAngles_y;
LShoulder_total(:,3) = LShoulderAngles_z;

LShoulder_total_max = max(LShoulder_total);
LShoulder_total_min = min(LShoulder_total);

for i = 1:length(LShoulder_total_max)
    LShoulder_total_ROM(1,i) = LShoulder_total_max(i)- LShoulder_total_min(i);
end

%ROM_Right_Shoulder = table([RShoulder_ph1_ROM(1);RShoulder_ph2_ROM(1);RShoulder_ph3_ROM(1);RShoulder_ph4_ROM(1)],[RShoulder_ph1_ROM(2);RShoulder_ph2_ROM(2);RShoulder_ph3_ROM(2);RShoulder_ph4_ROM(2)],[RShoulder_ph1_ROM(3);RShoulder_ph2_ROM(3);RShoulder_ph3_ROM(3);RShoulder_ph4_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Phase 1','Phase 2','Phase 3', 'Phase 4'})
ROM_Right_Shoulder = struct;
ROM_Right_Shoulder.Phase1_X=RShoulder_ph1_ROM(1);
ROM_Right_Shoulder.Phase1_Y=RShoulder_ph1_ROM(2);
ROM_Right_Shoulder.Phase1_Z=RShoulder_ph1_ROM(3);
ROM_Right_Shoulder.Phase2_X=RShoulder_ph2_ROM(1);
ROM_Right_Shoulder.Phase2_Y=RShoulder_ph2_ROM(2);
ROM_Right_Shoulder.Phase2_Z=RShoulder_ph2_ROM(3);
ROM_Right_Shoulder.Phase3_X=RShoulder_ph3_ROM(1);
ROM_Right_Shoulder.Phase3_Y=RShoulder_ph3_ROM(2);
ROM_Right_Shoulder.Phase3_Z=RShoulder_ph3_ROM(3);
ROM_Right_Shoulder.Phase4_X=RShoulder_ph4_ROM(1);
ROM_Right_Shoulder.Phase4_Y=RShoulder_ph4_ROM(2);
ROM_Right_Shoulder.Phase4_Z=RShoulder_ph4_ROM(3);

ROM_Right_Shoulder

%ROM_Left_Shoulder = table([LShoulder_ph1_ROM(1);LShoulder_ph2_ROM(1);LShoulder_ph3_ROM(1);LShoulder_ph4_ROM(1)],[LShoulder_ph1_ROM(2);LShoulder_ph2_ROM(2);LShoulder_ph3_ROM(2);LShoulder_ph4_ROM(2)],[LShoulder_ph1_ROM(3);LShoulder_ph2_ROM(3);LShoulder_ph3_ROM(3);LShoulder_ph4_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Phase 1','Phase 2','Phase 3', 'Phase 4'})
ROM_Left_Shoulder = struct;
ROM_Left_Shoulder.Phase1_X=LShoulder_ph1_ROM(1);
ROM_Left_Shoulder.Phase1_Y=LShoulder_ph1_ROM(2);
ROM_Left_Shoulder.Phase1_Z=LShoulder_ph1_ROM(3);
ROM_Left_Shoulder.Phase2_X=LShoulder_ph2_ROM(1);
ROM_Left_Shoulder.Phase2_Y=LShoulder_ph2_ROM(2);
ROM_Left_Shoulder.Phase2_Z=LShoulder_ph2_ROM(3);
ROM_Left_Shoulder.Phase3_X=LShoulder_ph3_ROM(1);
ROM_Left_Shoulder.Phase3_Y=LShoulder_ph3_ROM(2);
ROM_Left_Shoulder.Phase3_Z=LShoulder_ph3_ROM(3);
ROM_Left_Shoulder.Phase4_X=LShoulder_ph4_ROM(1);
ROM_Left_Shoulder.Phase4_Y=LShoulder_ph4_ROM(2);
ROM_Left_Shoulder.Phase4_Z=LShoulder_ph4_ROM(3);

ROM_Left_Shoulder

%TotalROM_Shoulder = table([RShoulder_total_ROM(1); LShoulder_total_ROM(1)],[RShoulder_total_ROM(2); LShoulder_total_ROM(2)],[RShoulder_total_ROM(3); LShoulder_total_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Right','Left'})
TotalROM_Shoulder = struct;
TotalROM_Shoulder.Right_X=RShoulder_total_ROM(1);
TotalROM_Shoulder.Right_Y=RShoulder_total_ROM(2);
TotalROM_Shoulder.Right_Z=RShoulder_total_ROM(3);
TotalROM_Shoulder.Left_X=LShoulder_total_ROM(1);
TotalROM_Shoulder.Left_Y=LShoulder_total_ROM(2);
TotalROM_Shoulder.Left_Z=LShoulder_total_ROM(3);

TotalROM_Shoulder
