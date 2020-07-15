%% Range of movement of Shoulder in Lateral Box Transfer

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

load('..\tests\data\input\dinamica08_B.mat')
             
Ts = 1/Fs;
t_total = (double(frames)*Ts);
t_1000 = 0:(t_total/1000):(t_total-(t_total/1000));

% Normal Ranges of Shoulder Motion
for i =1:(length(t_1000))
    flex_shoulder(1,i) = 130; % flexion
    ext_shoulder(1,i) = -30;   % extension
    abd_shoulder(1,i) = -50;   % abduction
    add_shoulder(1,i) = 30;   % adduction
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

LHEE_z_z = zscore(LHEE_z);
RHEE_z_z = zscore(RHEE_z);


t= 0:(length(LHEE)-1);
%plot(t, zscore(LHEE_z));

plot(t,LHEE_z_z);
hold on
%plot(t,zscore(RHEE_z));
plot(t,RHEE_z_z);

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
legend('Left Heel Z', 'Right Heel Z');
xlabel('Time in frames');
ylabel('Milimetres');
title('Segmentation Without Exoskeleton Trials')

%% Import model data for Shoulder Relative Angles

% Shoulder  
        RShoulderAngles_x = ModelData.Raw.(ModelOutput{9})(1,:); 
        RShoulderAngles_y = ModelData.Raw.(ModelOutput{9})(2,:);
        RShoulderAngles_z = ModelData.Raw.(ModelOutput{9})(3,:);
        

        LShoulderAngles_x = ModelData.Raw.(ModelOutput{10})(1,:);
        LShoulderAngles_y = ModelData.Raw.(ModelOutput{10})(2,:);
        LShoulderAngles_z = ModelData.Raw.(ModelOutput{10})(3,:);
        

%% Shoulder signals visualization before going through segmentation

figure(2)
frames = length(LHEE_z);
% (1:(frames/1000):frames)  si frames>1000
% (1:(frames/1000):frames+(frames/1000)) si frames<1000
% +(frames/1000)
plot(t_1000,RShoulderAngles_x(1:(frames/1000):frames)); 
hold on 
plot(t_1000, flex_shoulder,'--',t_1000, ext_shoulder,'--'); 
legend({'x', 'Normal range of flexion-extension'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
%ylim([-100 160]);
title('Right Shoulder angles: Flexo-Extension');

figure(3)
plot(t_1000, RShoulderAngles_y(1:(frames/1000):frames));
hold on 
plot(t_1000, abd_shoulder,'--',t_1000, add_shoulder,'--'); 
legend({'y', 'Normal range of abduction-adduction'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
ylim([-60 100]);
title('Right Shoulder angles: Abduction-Adduction');

figure(4)
plot(t_1000, RShoulderAngles_z(1:(frames/1000):frames));
legend({'z'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
%ylim([-70 50]);
title('Right Shoulder angles: Internal-External Rotation');

figure(5)
% (1:(frames/1000):frames)  si frames>1000
% (1:(frames/1000):frames+(frames/1000)) si frames<1000
% +(frames/1000)
plot(t_1000,LShoulderAngles_x(1:(frames/1000):frames)); 
hold on 
plot(t_1000, flex_shoulder,'--',t_1000, ext_shoulder,'--'); 
legend({'x', 'Normal range of flexion-extension'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
%ylim([-100 160]);
title('Left Shoulder angles: Flexo-Extension');

figure(6)
plot(t_1000, LShoulderAngles_y(1:(frames/1000):frames));
hold on 
plot(t_1000, abd_shoulder,'--',t_1000, add_shoulder,'--'); 
legend({'y', 'Normal range of abduction-adduction'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
ylim([-60 100]);
title('Left Shoulder angles: Abduction-Adduction');

figure(7)
plot(t_1000, LShoulderAngles_z(1:(frames/1000):frames));
legend({'z', 'Normal range of internal-external rotation'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
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

phase1_RShoulder_x = RShoulderAngles_x(idxO1:idxF2);  
phase1_RShoulder_y = RShoulderAngles_y(idxO1:idxF2);
phase1_RShoulder_z = RShoulderAngles_z(idxO1:idxF2);

phase2_RShoulder_x = RShoulderAngles_x(idxF2:idxF4); 
phase2_RShoulder_y = RShoulderAngles_y(idxF2:idxF4);
phase2_RShoulder_z = RShoulderAngles_z(idxF2:idxF4);

phase3_RShoulder_x = RShoulderAngles_x(idxF4:idxF5); 
phase3_RShoulder_y = RShoulderAngles_y(idxF4:idxF5);
phase3_RShoulder_z = RShoulderAngles_z(idxF4:idxF5);

% Left Shoulder

phase1_LShoulder_x = LShoulderAngles_x(idxO1:idxF2);
phase1_LShoulder_y = LShoulderAngles_y(idxO1:idxF2);
phase1_LShoulder_z = LShoulderAngles_z(idxO1:idxF2);

phase2_LShoulder_x = LShoulderAngles_x(idxF2:idxF4);
phase2_LShoulder_y = LShoulderAngles_y(idxF2:idxF4);
phase2_LShoulder_z = LShoulderAngles_z(idxF2:idxF4);

phase3_LShoulder_x = LShoulderAngles_x(idxF4:idxF5);
phase3_LShoulder_y = LShoulderAngles_y(idxF4:idxF5);
phase3_LShoulder_z = LShoulderAngles_z(idxF4:idxF5);

%% time vectors definition

t2 = 0:1:(length(phase3_RShoulder_x)-1);
t3 = 0:1:(length(phase1_RShoulder_x)-1);
t4 = 0:1:(length(phase2_RShoulder_x)-1);

% Right Shoulder figures: MOVEMENT SEGMENTED INTO THE 4 PHASES

figure(8)
subplot(3,1,1)
plot (t3, phase1_RShoulder_x, t3, phase1_RShoulder_y, '--',t3,phase1_RShoulder_z, '.');
legend('x','y','z')
xlabel('Time in frames'); 
ylabel('Degrees');
title('RIGHT ANGLES: Phase 1'); % Box Carrying from Sagittal to Frontal
subplot(3,1,2)
plot(t4, phase2_RShoulder_x,t4,phase2_RShoulder_y, '--',t4,phase2_RShoulder_z, '.');
legend('x','y','z')
xlabel('Time in frames'); 
ylabel('Degrees');
title('Phase 2'); % Box Carrying from Frontal to Sagittal
subplot(3,1,3)
plot(t2, phase3_RShoulder_x, t2, phase3_RShoulder_y,'--', t2, phase3_RShoulder_z, '.');
legend('x','y','z')
xlabel('Time in frames'); 
ylabel('Degrees');
title('Phase 3'); % Placing box sagittal plane again


% Left Shoulder figures: MOVEMENT SEGMENTED INTO THE 4 PHASES

figure(9)
subplot(3,1,1)
plot (t3, phase1_LShoulder_x, t3, phase1_LShoulder_y, '--',t3,phase1_LShoulder_z, '.');
legend('x','y','z')
xlabel('Time in frames'); 
ylabel('Degrees');
title('LEFT ANGLES: Phase 1'); % Box Carrying from Sagittal to Frontal
subplot(3,1,2)
plot(t4, phase2_LShoulder_x,t4, phase2_LShoulder_y, '--',t4, phase2_LShoulder_z, '.');
legend('x','y','z')
xlabel('Time in frames'); 
ylabel('Degrees');
title('Phase 2'); % Box Carrying from Frontal to Sagittal
subplot(3,1,3)
plot(t2, phase3_LShoulder_x, t2, phase3_LShoulder_y,'--', t2, phase3_LShoulder_z, '.');
legend('x','y','z')
xlabel('Time in frames'); 
ylabel('Degrees');
title('Phase 3'); % Placing box sagittal plane again


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


RShoulder_ph1_max = max(RShoulder_ph1);
RShoulder_ph1_min = min(RShoulder_ph1);

RShoulder_ph2_max = max(RShoulder_ph2);
RShoulder_ph2_min = min(RShoulder_ph2);

RShoulder_ph3_max = max(RShoulder_ph3);
RShoulder_ph3_min = min(RShoulder_ph3);


for i = 1:length(RShoulder_ph1_max)
    RShoulder_ph1_ROM(1,i) = RShoulder_ph1_max(i)- RShoulder_ph1_min(i);
end

for i = 1:length(RShoulder_ph1_max)
    RShoulder_ph2_ROM(1,i) = RShoulder_ph2_max(i)- RShoulder_ph2_min(i);
end

for i = 1:length(RShoulder_ph1_max)
    RShoulder_ph3_ROM(1,i) = RShoulder_ph3_max(i)- RShoulder_ph3_min(i);
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


LShoulder_ph1_max = max(LShoulder_ph1);
LShoulder_ph1_min = min(LShoulder_ph1);

LShoulder_ph2_max = max(LShoulder_ph2);
LShoulder_ph2_min = min(LShoulder_ph2);

LShoulder_ph3_max = max(LShoulder_ph3);
LShoulder_ph3_min = min(LShoulder_ph3);


for i = 1:length(RShoulder_ph1_max)
    LShoulder_ph1_ROM(1,i) = LShoulder_ph1_max(i)- LShoulder_ph1_min(i);
end

for i = 1:length(RShoulder_ph1_max)
    LShoulder_ph2_ROM(1,i) = LShoulder_ph2_max(i)- LShoulder_ph2_min(i);
end

for i = 1:length(RShoulder_ph1_max)
    LShoulder_ph3_ROM(1,i) = LShoulder_ph3_max(i)- LShoulder_ph3_min(i);
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

%ROM_Right_Shoulder = table([RShoulder_ph1_ROM(1);RShoulder_ph2_ROM(1);RShoulder_ph3_ROM(1)],[RShoulder_ph1_ROM(2);RShoulder_ph2_ROM(2);RShoulder_ph3_ROM(2)],[RShoulder_ph1_ROM(3);RShoulder_ph2_ROM(3);RShoulder_ph3_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Phase 1','Phase 2','Phase 3'})
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

ROM_Right_Shoulder

%ROM_Left_Shoulder = table([LShoulder_ph1_ROM(1);LShoulder_ph2_ROM(1);LShoulder_ph3_ROM(1)],[LShoulder_ph1_ROM(2);LShoulder_ph2_ROM(2);LShoulder_ph3_ROM(2)],[LShoulder_ph1_ROM(3);LShoulder_ph2_ROM(3);LShoulder_ph3_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Phase 1','Phase 2','Phase 3'})
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

