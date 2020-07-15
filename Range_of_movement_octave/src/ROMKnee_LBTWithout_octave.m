%% Range of movement of Knee in Lateral Box Transfer

% Code AUTHOR: Yaiza Benito Molpeceres. DATE: January-May 2020.
% Adapted to Octave by Guillermo Asín Prieto

% Without exoskeleton trials segmentation: din 42, 43, 47, 48, 49, 08, 21.
% Three phases shown in figure 1:
    % 1: Subject picks up the box in the sagittal plane and takes a step to rotate to the frontal plane (idxO1:idxF2)
    % 2: Subject deposits the box in the frontal plane and takes a step to rotate back to the sagittal plane(idxF2:idxF4)
    % 3: Subject deposits the box in the sagittal plane (idxF4:idxF5)

pkg load signal

clear all % Clear variables
close all % Close figures
clc

load('..\tests\data\input\dinamica42_B.mat')
             
Ts = 1/Fs;
t_total = (double(frames)*Ts);
t_1000 = 0:(t_total/1000):(t_total-(t_total/1000));

% Normal Ranges of Knee Motion
for i =1:(length(t_1000))
    flex_knee(1,i) = 130; % flexion
    ext_knee(1,i) = -15;   % extension
    intRot_hip(1,i) = 10;   % internal rotation
end

% Feet markers for signal segmentation
%LHEE_z = LTOE(:,3)'; 
%RHEE_z = RTOE(:,3)';   
LHEE_z = LHEE(:,3)';  
RHEE_z = RHEE(:,3)';    
  
%% Segmentation using LHEE, RHEE signals

% First turn

%[pks, locs] = findpeaks(LHEE_z,'minPeakProminence',10);
[pks, locs] = findpeaks(LHEE_z,'DoubleSided');   
%locs=locs(pks>80);                 
%pks=pks(pks>80);
for count=1:1:length(pks)
    if count==1
        signal_aux_L=LHEE_z(1:locs(count));
        signal_aux_R=LHEE_z(locs(count):locs(count+1));
    else if count==length(pks)
            signal_aux_L=LHEE_z(locs(count-1):locs(count));
            signal_aux_R=LHEE_z(locs(count):end);
        else
            signal_aux_L=LHEE_z(locs(count-1):locs(count));
            signal_aux_R=LHEE_z(locs(count):locs(count+1));
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
pks=pks(prominence>10);    
locs=locs(prominence>10);



%TF1 = islocalmin(LHEE_z, 'FlatSelection', 'first');
[TF1,locs_TF1] = findpeaks(1.01*max(LHEE_z(1:locs(1)))-LHEE_z(1:locs(1)),'DoubleSided');
TF1=1.01*max(LHEE_z(1:locs(1)))-TF1;

TF1_aux=zeros(1,length(TF1));
for count=1:1:length(TF1)-1
  if TF1(count+1)==TF1(count)
    TF1_aux(count+1)=1;
  end
end
TF1(find(TF1_aux))=[];
locs_TF1(find(TF1_aux))=[]; 


for count_aux=1:1:length(pks)-1
  %TF1 = islocalmin(LHEE_z, 'FlatSelection','first');
  [TF1_2,locs_TF1_2] = findpeaks(1.01*max(LHEE_z(locs(count_aux):locs(count_aux+1)))-LHEE_z(locs(count_aux):locs(count_aux+1)),'DoubleSided');
  TF1_2=1.01*max(LHEE_z(locs(count_aux):locs(count_aux+1)))-TF1_2;

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


%TF1 = islocalmin(LHEE_z, 'FlatSelection','first');
[TF1_2,locs_TF1_2] = findpeaks(1.01*max(LHEE_z(locs(end):end))-LHEE_z(locs(end):end),'DoubleSided');
TF1_2=1.01*max(LHEE_z(locs(end):end))-TF1_2;

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




%idx = find(TF1);
idx = locs_TF1;


flat2 = idx > locs(1);
idx_flat2 = find(flat2);

%[pks2, locs2] = findpeaks(RHEE_z,'minPeakProminence',10);
[pks2, locs2] = findpeaks(RHEE_z,'DoubleSided');   
%locs2=locs2(pks2>80);                 
%pks2=pks2(pks2>80);
for count=1:1:length(pks2)
    if count==1
        signal_aux_L=RHEE_z(1:locs2(count));
        signal_aux_R=RHEE_z(locs2(count):locs2(count+1));
    else if count==length(pks2)
            signal_aux_L=RHEE_z(locs2(count-1):locs2(count));
            signal_aux_R=RHEE_z(locs2(count):end);
        else
            signal_aux_L=RHEE_z(locs2(count-1):locs2(count));
            signal_aux_R=RHEE_z(locs2(count):locs2(count+1));
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
    prominence(count)=pks2(count)-real_min;
end
pks2=pks2(prominence>10);    
locs2=locs2(prominence>10);





%TF2 = islocalmin(RHEE_z, 'FlatSelection', 'first');
[TF2,locs_TF2] = findpeaks(1.01*max(RHEE_z(1:locs(1)))-RHEE_z(1:locs(1)),'DoubleSided');
TF2=1.01*max(RHEE_z(1:locs(1)))-TF2;

TF2_aux=zeros(1,length(TF2));
for count=1:1:length(TF2)-1
  if TF2(count+1)==TF2(count)
    TF2_aux(count+1)=1;
  end
end
TF2(find(TF2_aux))=[];
locs_TF2(find(TF2_aux))=[]; 


for count_aux=1:1:length(pks)-1
  %TF2 = islocalmin(RHEE_z, 'FlatSelection','first');
  [TF2_2,locs_TF2_2] = findpeaks(1.01*max(RHEE_z(locs(count_aux):locs(count_aux+1)))-RHEE_z(locs(count_aux):locs(count_aux+1)),'DoubleSided');
  TF2_2=1.01*max(RHEE_z(locs(count_aux):locs(count_aux+1)))-TF2_2;

  TF2_aux_2=zeros(1,length(TF2_2));
  for count=1:1:length(TF2_2)-1
    if TF2_2(count+1)==TF2_2(count)
      TF2_aux_2(count+1)=1;
    end
  end
  TF2_2(find(TF2_aux_2))=[];
  locs_TF2_2(find(TF2_aux_2))=[];
  TF2 = [TF2,TF2_2];
  locs_TF2 = [locs_TF2,locs_TF2_2+locs(count_aux)];
endfor


%TF2 = islocalmin(RHEE_z, 'FlatSelection','first');
[TF2_2,locs_TF2_2] = findpeaks(1.01*max(RHEE_z(locs(end):end))-RHEE_z(locs(end):end),'DoubleSided');
TF2_2=1.01*max(RHEE_z(locs(end):end))-TF2_2;

TF2_aux_2=zeros(1,length(TF2_2));
for count=1:1:length(TF2_2)-1
  if TF2_2(count+1)==TF2_2(count)
    TF2_aux_2(count+1)=1;
  end
end
TF2_2(find(TF2_aux_2))=[];
locs_TF2_2(find(TF2_aux_2))=[]; 

TF2 = [TF2,TF2_2];
locs_TF2 = [locs_TF2,locs_TF2_2+locs(end)];

for count=1:1:length(TF2)-1
    for count_aux=1:1:length(pks)
        if TF2(count)==pks(count_aux)
          TF2_aux(count)=1;
        end
    end
end

TF2(find(TF2_aux))=[];
locs_TF2(find(TF2_aux))=[];









%idx2 = find(TF2);
idx2 = locs_TF2;


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
t= (0:(length(LHEE)-1))*Ts;

RHEE_z_z = zscore(RHEE_z);
LHEE_z_z = zscore(LHEE_z);

plot(t,LHEE_z_z);
hold on
plot(t,RHEE_z_z);

%plot(t, zscore(LHEE_z));
%hold on
%plot(t,zscore(RHEE_z));

%plot(t,zscore(LHEE_z),'r*', 'MarkerIndices', idxO1);
plot(t(idxO1),LHEE_z_z(idxO1),'r*');

%plot(t,zscore(LHEE_z),'r*', 'MarkerIndices', idxF2); % 42 idxF2
plot(t(idxF2),LHEE_z_z(idxF2),'r*');

%plot(t,zscore(LHEE_z),'r*', 'MarkerIndices', idxF5);
plot(t(idxF5),LHEE_z_z(idxF5),'r*');

%plot(t,zscore(RHEE_z),'r*', 'MarkerIndices', idxO1);
plot(t(idxO1),RHEE_z_z(idxO1),'r*');

%plot(t,zscore(RHEE_z),'r*', 'MarkerIndices', idxF4); % 42 idxF4
plot(t(idxF4),RHEE_z_z(idxF4),'r*');

%plot(t,zscore(RHEE_z),'r*', 'MarkerIndices', idxF5);
plot(t(idxF5),RHEE_z_z(idxF5),'r*');

hold off
xlabel('Time in seconds'),
ylabel('Angular displacement in degrees');
legend('Left Heel Z', 'Right Heel Z');
title('Segmentation Without Exoskeleton Trials')


%% Import model data for Knee Relative Angles

% Knee  
        RKneeAngles_x = ModelData.Raw.(ModelOutput{3})(1,:); 
        RKneeAngles_y = ModelData.Raw.(ModelOutput{3})(2,:);
        RKneeAngles_z = ModelData.Raw.(ModelOutput{3})(3,:);
        

        LKneeAngles_x = ModelData.Raw.(ModelOutput{4})(1,:);
        LKneeAngles_y = ModelData.Raw.(ModelOutput{4})(2,:);
        LKneeAngles_z = ModelData.Raw.(ModelOutput{4})(3,:);
        

%% Knee signals visualization before going through segmentation

figure(2)
frames = length(LHEE_z);
% (1:(frames/1000):frames)  si frames>1000
% (1:(frames/1000):frames+(frames/1000)) si frames<1000
% +(frames/1000)
plot(t_1000,RKneeAngles_x(1:(frames/1000):frames)); 
hold on 
plot(t_1000, flex_knee,'--',t_1000, ext_knee,'--'); 
legend({'x', 'Normal range of flexion-extension'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
%ylim([-100 160]);
title('Right Knee angles: Flexo-Extension');

figure(3)
plot(t_1000, RKneeAngles_y(1:(frames/1000):frames));
legend({'y'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');  
%ylim([-70 50]);
title('Right Knee angles: Abduction-Adduction');

figure(4)
plot(t_1000, RKneeAngles_z(1:(frames/1000):frames));
hold on 
plot(t_1000, intRot_hip,'--'); 
legend({'z', 'Normal range of internal rotation'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil  
ylabel('Degrees');
ylim([-70 20]);
title('Right Knee angles: Internal-External Rotation');

figure(5)
% (1:(frames/1000):frames)  si frames>1000
% (1:(frames/1000):frames+(frames/1000)) si frames<1000
% +(frames/1000)
plot(t_1000,LKneeAngles_x(1:(frames/1000):frames)); 
hold on 
plot(t_1000, flex_knee,'--',t_1000, ext_knee,'--'); 
legend({'x', 'Normal range of flexion-extension'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
%ylim([-100 160]);
title('Left Knee angles: Flexo-Extension');

figure(6)
plot(t_1000, LKneeAngles_y(1:(frames/1000):frames));
legend({'y'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
%ylim([-70 50]);
title('Left Knee angles: Abduction-Adduction');

figure(7)
plot(t_1000, LKneeAngles_z(1:(frames/1000):frames));
hold on 
plot(t_1000, intRot_hip,'--'); 
legend({'z', 'Normal range of internal rotation'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
ylim([-70 20]);
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

phase2_RKnee_x = RKneeAngles_x(idxF1:idxF3); 
phase2_RKnee_y = RKneeAngles_y(idxF1:idxF3);
phase2_RKnee_z = RKneeAngles_z(idxF1:idxF3);

phase3_RKnee_x = RKneeAngles_x(idxF3:idxF5); 
phase3_RKnee_y = RKneeAngles_y(idxF3:idxF5);
phase3_RKnee_z = RKneeAngles_z(idxF3:idxF5);

% Left Knee

phase1_LKnee_x = LKneeAngles_x(idxO1:idxF1);
phase1_LKnee_y = LKneeAngles_y(idxO1:idxF1);
phase1_LKnee_z = LKneeAngles_z(idxO1:idxF1);

phase2_LKnee_x = LKneeAngles_x(idxF1:idxF3);
phase2_LKnee_y = LKneeAngles_y(idxF1:idxF3);
phase2_LKnee_z = LKneeAngles_z(idxF1:idxF3);

phase3_LKnee_x = LKneeAngles_x(idxF3:idxF5);
phase3_LKnee_y = LKneeAngles_y(idxF3:idxF5);
phase3_LKnee_z = LKneeAngles_z(idxF3:idxF5);

%% time vectors definition

t2 = 0:1:(length(phase3_RKnee_x)-1);
t3 = 0:1:(length(phase1_RKnee_x)-1);
t4 = 0:1:(length(phase2_RKnee_x)-1);

% Right Knee figures: MOVEMENT SEGMENTED INTO THE 4 PHASES

figure(8)
subplot(3,1,1)
plot (t3, phase1_RKnee_x, t3, phase1_RKnee_y, '--',t3,phase1_RKnee_z, '.');
legend('x','y','z')
xlabel('Time in frames');
ylabel('Degrees');
title('RIGHT ANGLES: Phase 1'); % Box Carrying from Sagittal to Frontal
subplot(3,1,2)
plot(t4, phase2_RKnee_x,t4,phase2_RKnee_y, '--',t4,phase2_RKnee_z, '.');
legend('x','y','z')
xlabel('Time in frames');
ylabel('Degrees');
title('Phase 2'); % Box Carrying from Frontal to Sagittal
subplot(3,1,3)
plot(t2, phase3_RKnee_x, t2, phase3_RKnee_y,'--', t2, phase3_RKnee_z, '.');
legend('x','y','z')
xlabel('Time in frames');
ylabel('Degrees');
title('Phase 3'); % Placing box sagittal plane again


% Left Knee figures: MOVEMENT SEGMENTED INTO THE 4 PHASES

figure(9)
subplot(3,1,1)
plot (t3, phase1_LKnee_x, t3, phase1_LKnee_y, '--',t3,phase1_LKnee_z, '.');
legend('x','y','z')
xlabel('Time in frames');
ylabel('Degrees');
title('LEFT ANGLES: Phase 1'); % Box Carrying from Sagittal to Frontal
subplot(3,1,2)
plot(t4, phase2_LKnee_x,t4, phase2_LKnee_y, '--',t4, phase2_LKnee_z, '.');
legend('x','y','z')
xlabel('Time in frames');
ylabel('Degrees');
title('Phase 2'); % Box Carrying from Frontal to Sagittal
subplot(3,1,3)
plot(t2, phase3_LKnee_x, t2, phase3_LKnee_y,'--', t2, phase3_LKnee_z, '.');
legend('x','y','z')
xlabel('Time in frames');
ylabel('Degrees');
title('Phase 3'); % Placing box sagittal plane again


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


RKnee_ph1_max = max(RKnee_ph1);
RKnee_ph1_min = min(RKnee_ph1);

RKnee_ph2_max = max(RKnee_ph2);
RKnee_ph2_min = min(RKnee_ph2);

RKnee_ph3_max = max(RKnee_ph3);
RKnee_ph3_min = min(RKnee_ph3);


for i = 1:length(RKnee_ph1_max)
    RKnee_ph1_ROM(1,i) = RKnee_ph1_max(i)- RKnee_ph1_min(i);
end

for i = 1:length(RKnee_ph1_max)
    RKnee_ph2_ROM(1,i) = RKnee_ph2_max(i)- RKnee_ph2_min(i);
end

for i = 1:length(RKnee_ph1_max)
    RKnee_ph3_ROM(1,i) = RKnee_ph3_max(i)- RKnee_ph3_min(i);
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


LKnee_ph1_max = max(LKnee_ph1);
LKnee_ph1_min = min(LKnee_ph1);

LKnee_ph2_max = max(LKnee_ph2);
LKnee_ph2_min = min(LKnee_ph2);

LKnee_ph3_max = max(LKnee_ph3);
LKnee_ph3_min = min(LKnee_ph3);


for i = 1:length(RKnee_ph1_max)
    LKnee_ph1_ROM(1,i) = LKnee_ph1_max(i)- LKnee_ph1_min(i);
end

for i = 1:length(RKnee_ph1_max)
    LKnee_ph2_ROM(1,i) = LKnee_ph2_max(i)- LKnee_ph2_min(i);
end

for i = 1:length(RKnee_ph1_max)
    LKnee_ph3_ROM(1,i) = LKnee_ph3_max(i)- LKnee_ph3_min(i);
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

%ROM_Right_Knee = table([RKnee_ph1_ROM(1);RKnee_ph2_ROM(1);RKnee_ph3_ROM(1)],[RKnee_ph1_ROM(2);RKnee_ph2_ROM(2);RKnee_ph3_ROM(2)],[RKnee_ph1_ROM(3);RKnee_ph2_ROM(3);RKnee_ph3_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Phase 1','Phase 2','Phase 3'})
ROM_Right_Knee = struct;
ROM_Right_Knee.Phase1_X=RKnee_ph1_ROM(1);
ROM_Right_Knee.Phase1_Y=RKnee_ph1_ROM(2);
ROM_Right_Knee.Phase1_Z=RKnee_ph1_ROM(3);
ROM_Right_Knee.Phase2_X=RKnee_ph2_ROM(1);
ROM_Right_Knee.Phase2_Y=RKnee_ph2_ROM(2);
ROM_Right_Knee.Phase2_Z=RKnee_ph2_ROM(3);
ROM_Right_Knee.Phase3_X=RKnee_ph3_ROM(1);
ROM_Right_Knee.Phase3_Y=RKnee_ph3_ROM(2);
ROM_Right_Knee.Phase3_Z=RKnee_ph3_ROM(3);

ROM_Right_Knee

%ROM_Left_Knee = table([LKnee_ph1_ROM(1);LKnee_ph2_ROM(1);LKnee_ph3_ROM(1)],[LKnee_ph1_ROM(2);LKnee_ph2_ROM(2);LKnee_ph3_ROM(2)],[LKnee_ph1_ROM(3);LKnee_ph2_ROM(3);LKnee_ph3_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Phase 1','Phase 2','Phase 3'})
ROM_Left_Knee = struct;
ROM_Left_Knee.Phase1_X=LKnee_ph1_ROM(1);
ROM_Left_Knee.Phase1_Y=LKnee_ph1_ROM(2);
ROM_Left_Knee.Phase1_Z=LKnee_ph1_ROM(3);
ROM_Left_Knee.Phase2_X=LKnee_ph2_ROM(1);
ROM_Left_Knee.Phase2_Y=LKnee_ph2_ROM(2);
ROM_Left_Knee.Phase2_Z=LKnee_ph2_ROM(3);
ROM_Left_Knee.Phase3_X=LKnee_ph3_ROM(1);
ROM_Left_Knee.Phase3_Y=LKnee_ph3_ROM(2);
ROM_Left_Knee.Phase3_Z=LKnee_ph3_ROM(3);

ROM_Left_Knee

%TotalROM_Knee = table([RKnee_total_ROM(1); LKnee_total_ROM(1)],[RKnee_total_ROM(2); LKnee_total_ROM(2)],[RKnee_total_ROM(3); LKnee_total_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Right','Left'})
TotalROM_Knee = struct;
TotalROM_Knee.Right_X=RKnee_total_ROM(1);
TotalROM_Knee.Right_Y=RKnee_total_ROM(2);
TotalROM_Knee.Right_Z=RKnee_total_ROM(3);
TotalROM_Knee.Left_X=LKnee_total_ROM(1);
TotalROM_Knee.Left_Y=LKnee_total_ROM(2);
TotalROM_Knee.Left_Z=LKnee_total_ROM(3);

TotalROM_Knee
