%% Range of movement of Ankle in Lateral Box Transfer

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

% Normal Ranges of Ankle Motion
for i =1:(length(t_1000))
    flex_ankle(1,i) = 45; % flexion
    ext_ankle(1,i) = -20;   % extension
    supination_ankle(1,i) = 20;  % supination
    pronation_ankle(1,i) = -20;  % pronation
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
t= 0:(length(LHEE)-1);

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
plot(t(idxF4),RHEE_z_z(idxF4),'r*');

hold off
legend('Left Heel Z', 'Right Heel Z');
title('Segmentation Without Exoskeleton Trials')

%% Import model data for Ankle Relative Angles

% Ankle  
        RAnkleAngles_x = ModelData.Raw.(ModelOutput{5})(1,:); 
        RAnkleAngles_y = ModelData.Raw.(ModelOutput{5})(2,:);
        RAnkleAngles_z = ModelData.Raw.(ModelOutput{5})(3,:);
        

        LAnkleAngles_x = ModelData.Raw.(ModelOutput{6})(1,:);
        LAnkleAngles_y = ModelData.Raw.(ModelOutput{6})(2,:);
        LAnkleAngles_z = ModelData.Raw.(ModelOutput{6})(3,:);
        

%% Ankle signals visualization before going through segmentation

figure(2)
frames = length(LHEE_z);
% (1:(frames/1000):frames)  si frames>1000
% (1:(frames/1000):frames+(frames/1000)) si frames<1000
% +(frames/1000)
plot(t_1000,RAnkleAngles_x(1:(frames/1000):frames+(frames/1000))); 
hold on 
plot(t_1000, flex_ankle,'--',t_1000, ext_ankle,'--'); 
legend({'x', 'Normal range of flexion-extension'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
ylim([-30 50]);
title('Right Ankle angles: Flexo-Extension');

figure(3)
plot(t_1000, RAnkleAngles_y(1:(frames/1000):frames+(frames/1000)));
hold on 
plot(t_1000, supination_ankle,'--',t_1000, pronation_ankle,'--'); % ver cual es la positiva
legend({'y', 'Normal range of supination-pronation'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
ylim([-30 30]);
title('Right Ankle angles: Supination-Pronation');

figure(4)
plot(t_1000, RAnkleAngles_z(1:(frames/1000):frames+(frames/1000)));
legend({'z'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
%ylim([-70 50]);
title('Right Ankle angles: Internal-External Rotation');

figure(5)
% (1:(frames/1000):frames)  si frames>1000
% (1:(frames/1000):frames+(frames/1000)) si frames<1000
% +(frames/1000)
plot(t_1000,LAnkleAngles_x(1:(frames/1000):frames+(frames/1000))); 
hold on 
plot(t_1000, flex_ankle,'--',t_1000, ext_ankle,'--'); 
legend({'x', 'Normal range of flexion-extension'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
ylim([-30 50]);;
title('Left Ankle angles: Flexo-Extension');

figure(6)
plot(t_1000, LAnkleAngles_y(1:(frames/1000):frames+(frames/1000)));
hold on 
plot(t_1000, supination_ankle,'--',t_1000, pronation_ankle,'--'); % ver cual es la positiva
legend({'y', 'Normal range of supination-pronation'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
ylim([-30 30]);
title('Left Ankle angles: Supination-Pronation');

figure(7)
plot(t_1000, LAnkleAngles_z(1:(frames/1000):frames+(frames/1000)));
legend({'z'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
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

phase2_RAnkle_x = RAnkleAngles_x(idxF1:idxF3); 
phase2_RAnkle_y = RAnkleAngles_y(idxF1:idxF3);
phase2_RAnkle_z = RAnkleAngles_z(idxF1:idxF3);

phase3_RAnkle_x = RAnkleAngles_x(idxF3:idxF5); 
phase3_RAnkle_y = RAnkleAngles_y(idxF3:idxF5);
phase3_RAnkle_z = RAnkleAngles_z(idxF3:idxF5);

% Left Ankle

phase1_LAnkle_x = LAnkleAngles_x(idxO1:idxF1);
phase1_LAnkle_y = LAnkleAngles_y(idxO1:idxF1);
phase1_LAnkle_z = LAnkleAngles_z(idxO1:idxF1);

phase2_LAnkle_x = LAnkleAngles_x(idxF1:idxF3);
phase2_LAnkle_y = LAnkleAngles_y(idxF1:idxF3);
phase2_LAnkle_z = LAnkleAngles_z(idxF1:idxF3);

phase3_LAnkle_x = LAnkleAngles_x(idxF3:idxF5);
phase3_LAnkle_y = LAnkleAngles_y(idxF3:idxF5);
phase3_LAnkle_z = LAnkleAngles_z(idxF3:idxF5);

%% time vectors definition

t2 = 0:1:(length(phase3_RAnkle_x)-1);
t3 = 0:1:(length(phase1_RAnkle_x)-1);
t4 = 0:1:(length(phase2_RAnkle_x)-1);

% Right Ankle figures: MOVEMENT SEGMENTED INTO THE 4 PHASES

figure(8)
subplot(3,1,1)
plot (t3, phase1_RAnkle_x, t3, phase1_RAnkle_y, '--',t3,phase1_RAnkle_z, '.');
legend('x','y','z')
xlabel('Time in frames'); 
ylabel('Degrees');
title('RIGHT ANGLES: Phase 1'); % Box Carrying from Sagittal to Frontal
subplot(3,1,2)
plot(t4, phase2_RAnkle_x,t4,phase2_RAnkle_y, '--',t4,phase2_RAnkle_z, '.');
legend('x','y','z')
xlabel('Time in frames'); 
ylabel('Degrees');
title('Phase 2'); % Box Carrying from Frontal to Sagittal
subplot(3,1,3)
plot(t2, phase3_RAnkle_x, t2, phase3_RAnkle_y,'--', t2, phase3_RAnkle_z, '.');
legend('x','y','z')
xlabel('Time in frames'); 
ylabel('Degrees');
title('Phase 3'); % Placing box sagittal plane again


% Left Ankle figures: MOVEMENT SEGMENTED INTO THE 4 PHASES

figure(9)
subplot(3,1,1)
plot (t3, phase1_LAnkle_x, t3, phase1_LAnkle_y, '--',t3,phase1_LAnkle_z, '.');
legend('x','y','z')
xlabel('Time in frames');  
ylabel('Degrees');
title('LEFT ANGLES: Phase 1'); % Box Carrying from Sagittal to Frontal
subplot(3,1,2)
plot(t4, phase2_LAnkle_x,t4, phase2_LAnkle_y, '--',t4, phase2_LAnkle_z, '.');
legend('x','y','z')
xlabel('Time in frames'); 
ylabel('Degrees');
title('Phase 2'); % Box Carrying from Frontal to Sagittal
subplot(3,1,3)
plot(t2, phase3_LAnkle_x, t2, phase3_LAnkle_y,'--', t2, phase3_LAnkle_z, '.');
legend('x','y','z')
xlabel('Time in frames'); 
ylabel('Degrees');
title('Phase 3'); % Placing box sagittal plane again


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


RAnkle_ph1_max = max(RAnkle_ph1);
RAnkle_ph1_min = min(RAnkle_ph1);

RAnkle_ph2_max = max(RAnkle_ph2);
RAnkle_ph2_min = min(RAnkle_ph2);

RAnkle_ph3_max = max(RAnkle_ph3);
RAnkle_ph3_min = min(RAnkle_ph3);


for i = 1:length(RAnkle_ph1_max)
    RAnkle_ph1_ROM(1,i) = RAnkle_ph1_max(i)- RAnkle_ph1_min(i);
end

for i = 1:length(RAnkle_ph1_max)
    RAnkle_ph2_ROM(1,i) = RAnkle_ph2_max(i)- RAnkle_ph2_min(i);
end

for i = 1:length(RAnkle_ph1_max)
    RAnkle_ph3_ROM(1,i) = RAnkle_ph3_max(i)- RAnkle_ph3_min(i);
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


LAnkle_ph1_max = max(LAnkle_ph1);
LAnkle_ph1_min = min(LAnkle_ph1);

LAnkle_ph2_max = max(LAnkle_ph2);
LAnkle_ph2_min = min(LAnkle_ph2);

LAnkle_ph3_max = max(LAnkle_ph3);
LAnkle_ph3_min = min(LAnkle_ph3);


for i = 1:length(RAnkle_ph1_max)
    LAnkle_ph1_ROM(1,i) = LAnkle_ph1_max(i)- LAnkle_ph1_min(i);
end

for i = 1:length(RAnkle_ph1_max)
    LAnkle_ph2_ROM(1,i) = LAnkle_ph2_max(i)- LAnkle_ph2_min(i);
end

for i = 1:length(RAnkle_ph1_max)
    LAnkle_ph3_ROM(1,i) = LAnkle_ph3_max(i)- LAnkle_ph3_min(i);
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

ROM_Right_Ankle = table([RAnkle_ph1_ROM(1);RAnkle_ph2_ROM(1);RAnkle_ph3_ROM(1)],[RAnkle_ph1_ROM(2);RAnkle_ph2_ROM(2);RAnkle_ph3_ROM(2)],[RAnkle_ph1_ROM(3);RAnkle_ph2_ROM(3);RAnkle_ph3_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Phase 1','Phase 2','Phase 3'})
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

ROM_Right_Ankle

%ROM_Left_Ankle = table([LAnkle_ph1_ROM(1);LAnkle_ph2_ROM(1);LAnkle_ph3_ROM(1)],[LAnkle_ph1_ROM(2);LAnkle_ph2_ROM(2);LAnkle_ph3_ROM(2)],[LAnkle_ph1_ROM(3);LAnkle_ph2_ROM(3);LAnkle_ph3_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Phase 1','Phase 2','Phase 3'})
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

ROM_Left_Ankle

TotalROM_Ankle = table([RAnkle_total_ROM(1); LAnkle_total_ROM(1)],[RAnkle_total_ROM(2); LAnkle_total_ROM(2)],[RAnkle_total_ROM(3); LAnkle_total_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Right','Left'})

TotalROM_Ankle = struct;
TotalROM_Ankle.Right_X=RAnkle_total_ROM(1);
TotalROM_Ankle.Right_Y=RAnkle_total_ROM(2);
TotalROM_Ankle.Right_Z=RAnkle_total_ROM(3);
TotalROM_Ankle.Left_X=LAnkle_total_ROM(1);
TotalROM_Ankle.Left_Y=LAnkle_total_ROM(2);
TotalROM_Ankle.Left_Z=LAnkle_total_ROM(3);

TotalROM_Ankle