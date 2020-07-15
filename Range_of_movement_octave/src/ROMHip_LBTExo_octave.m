%% Range of movement of Hip in Lateral Box Transfer

% Code AUTHOR: Yaiza Benito Molpeceres. DATE: January-May 2020.
% Adapted to Octave by Guillermo Asín Prieto

% Exoskeleton trials segmentation: 54, 55, 56, 62.
% Five phases shown in figure 1:
    % 1: Subject picks up the box in the sagittal plane (idxO1:idxF1)
    % 2: Subject takes some steps to rotate to the frontal plane (idxF1:idxF2)
    % 3: Subject deposits the box in the frontal plane (idxF2:idxF3)
    % 4: Subject takes some steps to rotate back to the sagittal plane (idxF3:idxF4)
    % 5: Subject deposits the box in the sagittal plane (idxF4:idxF5)


pkg load signal

clear all % Clear variables
close all % Close figures
clc

load('..\tests\data\input\dinamica55_B.mat')
           
Ts = 1/Fs;
t_total = (double(frames)*Ts);
t_1000 = 0:(t_total/1000):(t_total-(t_total/1000));

% Normal Ranges of Hip Motion
for i =1:(length(t_1000))
    flex_hip(1,i) = 130; % flexion
    ext_hip(1,i) = -30;   % extension
    abd_hip(1,i) = -50;   % abduction
    add_hip(1,i) = 30;   % adduction
    intRot_hip(1,i) = 40;   % internal rotation
    extRot_hip(1,i) = -45;   % external rotation
end

% Feet markers for signal segmentation
LTOE_z = LTOE(:,3)'; 
RTOE_z = RTOE(:,3)';   
%LHEE_z = LHEE(:,3)';  
%RHEE_z = RHEE(:,3)'; 

t= 0:(length(LTOE)-1);
%% Segmentation using LTOE, RTOE signals

% Left Toe

%[pks, locs] = findpeaks(LTOE_z,'minPeakProminence',10); 
[pks, locs] = findpeaks(LTOE_z,'DoubleSided');   
%locs=locs(pks>80);                 
%pks=pks(pks>80);
for count=1:1:length(pks)
    if count==1
        signal_aux_L=LTOE_z(1:locs(count));
        signal_aux_R=LTOE_z(locs(count):locs(count+1));
    else if count==length(pks)
            signal_aux_L=LTOE_z(locs(count-1):locs(count));
            signal_aux_R=LTOE_z(locs(count):end);
        else
            signal_aux_L=LTOE_z(locs(count-1):locs(count));
            signal_aux_R=LTOE_z(locs(count):locs(count+1));
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



e = locs < 500 | locs > 800;
f = find(e);
locs = locs(f);
pks = pks(f);
%TF1 = islocalmin(LTOE_z, 'FlatSelection', 'first');
%idx = find(TF1);
%flat = idx < locs(1);
%idx_flat = find(flat);

%flat2 = idx < locs(4)& idx > locs(3);
%idx_flat2 = find(flat2);

%flat3 = idx > locs(6);
%idx_flat3 = find(flat3);

% Right Toe 

%[pks2, locs2] = findpeaks(RTOE_z,'minPeakProminence',10); 
[pks2, locs2] = findpeaks(RTOE_z,'DoubleSided');   
%locs2=locs2(pks2>80);                 
%pks2=pks2(pks2>80);
for count=1:1:length(pks2)
    if count==1
        signal_aux_L=RTOE_z(1:locs2(count));
        signal_aux_R=RTOE_z(locs2(count):locs2(count+1));
    else if count==length(pks2)
            signal_aux_L=RTOE_z(locs2(count-1):locs2(count));
            signal_aux_R=RTOE_z(locs2(count):end);
        else
            signal_aux_L=RTOE_z(locs2(count-1):locs2(count));
            signal_aux_R=RTOE_z(locs2(count):locs2(count+1));
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





%TF2 = islocalmin(RTOE_z, 'FlatSelection', 'first');
[TF2,locs_TF2] = findpeaks(1.01*max(RTOE_z(1:locs(1)))-RTOE_z(1:locs(1)),'DoubleSided');
TF2=1.01*max(RTOE_z(1:locs(1)))-TF2;

TF2_aux=zeros(1,length(TF2));
for count=1:1:length(TF2)-1
  if TF2(count+1)==TF2(count)
    TF2_aux(count+1)=1;
  end
end
TF2(find(TF2_aux))=[];
locs_TF2(find(TF2_aux))=[]; 


for count_aux=1:1:length(pks)-1
  %TF2 = islocalmin(RTOE_z, 'FlatSelection','first');
  [TF2_2,locs_TF2_2] = findpeaks(1.01*max(RTOE_z(locs(count_aux):locs(count_aux+1)))-RTOE_z(locs(count_aux):locs(count_aux+1)),'DoubleSided');
  TF2_2=1.01*max(RTOE_z(locs(count_aux):locs(count_aux+1)))-TF2_2;

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


%TF2 = islocalmin(RTOE_z, 'FlatSelection','first');
[TF2_2,locs_TF2_2] = findpeaks(1.01*max(RTOE_z(locs(end):end))-RTOE_z(locs(end):end),'DoubleSided');
TF2_2=1.01*max(RTOE_z(locs(end):end))-TF2_2;

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


flat4 = idx2 < locs2(1);
idx_flat4 = find(flat4);

flat5 = idx2 > locs2(length(locs2));
idx_flat5 = find(flat5);
flat6 = idx2 > locs2(3)& idx2 < locs2(4);
idx_flat6 = find(flat6);

% Indexes definiton for signal segmentation
idxO1 = 1;
idxF1 = idx2(idx_flat4(length(idx_flat4))); %idxO2
idxF2 = idx2(idx_flat6(1)); %idxO3
idxF3 = idx2(idx_flat6(length(idx_flat6))); %idxO4
idxF4 = idx2(idx_flat5(1)); %idxO5
idxF5 = length(LTOE_z);

% Visualization of the segmentation
figure(1)
t= 0:(length(LTOE)-1);
RTOE_z_z = zscore(RTOE_z);
LTOE_z_z = zscore(LTOE_z);

plot(t,LTOE_z_z);
hold on
plot(t,RTOE_z_z);

%plot(t, zscore(LTOE_z));
%hold on
%plot(t,zscore(RTOE_z));

%plot(t,zscore(LTOE_z),'r*', 'MarkerIndices', idxO1);
plot(t(idxO1),LTOE_z_z(idxO1),'r*');

%plot(t,zscore(RTOE_z),'r*', 'MarkerIndices', idxF2);
plot(t(idxF2),RTOE_z_z(idxF2),'r*');

%plot(t,zscore(LTOE_z),'r*', 'MarkerIndices', idxF3);
plot(t(idxF3),LTOE_z_z(idxF3),'r*');

%plot(t,zscore(LTOE_z),'r*', 'MarkerIndices', idxF5);
plot(t(idxF5),LTOE_z_z(idxF5),'r*');

%plot(t,zscore(RTOE_z),'r*', 'MarkerIndices', idxO1);
plot(t(idxO1),RTOE_z_z(idxO1),'r*');

%plot(t,zscore(RTOE_z),'r*', 'MarkerIndices', idxF1);
plot(t(idxF1),RTOE_z_z(idxF1),'r*');

%plot(t,zscore(LTOE_z),'r*', 'MarkerIndices', idxF4);
plot(t(idxF4),LTOE_z_z(idxF4),'r*');

%plot(t,zscore(RTOE_z),'r*', 'MarkerIndices', idxF5);
plot(t(idxF5),RTOE_z_z(idxF5),'r*');

hold off
legend('Left Toe Z', 'Right Toe Z');
title('Segmentation Exoskeleton Trials')


%% Import model data for Hip Relative Angles

% Hip  
        RHipAngles_x = ModelData.Raw.(ModelOutput{1})(1,:); 
        RHipAngles_y = ModelData.Raw.(ModelOutput{1})(2,:);
        RHipAngles_z = ModelData.Raw.(ModelOutput{1})(3,:);
        

        LHipAngles_x = ModelData.Raw.(ModelOutput{2})(1,:);
        LHipAngles_y = ModelData.Raw.(ModelOutput{2})(2,:);
        LHipAngles_z = ModelData.Raw.(ModelOutput{2})(3,:);
        

%% Hip signals visualization before going through segmentation

figure(4)
frames = length(LTOE_z);
% (1:(frames/1000):frames)  si frames>1000
% (1:(frames/1000):frames+(frames/1000)) si frames<1000
% +(frames/1000)
plot(t_1000,RHipAngles_x(1:(frames/1000):frames)); 
hold on 
plot(t_1000, flex_hip,'--',t_1000, ext_hip,'--'); 
legend({'x', 'Normal range of flexion-extension'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
%ylim([-100 160]);
title('Right Hip angles: Flexo-Extension');

figure(5)
plot(t_1000, RHipAngles_y(1:(frames/1000):frames));
hold on 
plot(t_1000, abd_hip,'--',t_1000, add_hip,'--'); 
legend({'y', 'Normal range of abduction-adduction'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
ylim([-60 40]);
title('Right Hip angles: Abduction-Adduction');

figure(6)
plot(t_1000, RHipAngles_z(1:(frames/1000):frames));
hold on 
plot(t_1000, intRot_hip,'--',t_1000, extRot_hip,'--'); 
legend({'z', 'Normal range of internal-external rotation'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
%ylim([-70 50]);
title('Right Hip angles: Internal-External Rotation');

figure(7)
% (1:(frames/1000):frames)  si frames>1000
% (1:(frames/1000):frames+(frames/1000)) si frames<1000
% +(frames/1000)
plot(t_1000,LHipAngles_x(1:(frames/1000):frames)); 
hold on 
plot(t_1000, flex_hip,'--',t_1000, ext_hip,'--'); 
legend({'x', 'Normal range of flexion-extension'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
%ylim([-100 160]);
title('Left Hip angles: Flexo-Extension');

figure(8)
plot(t_1000, LHipAngles_y(1:(frames/1000):frames));
hold on 
plot(t_1000, abd_hip,'--',t_1000, add_hip,'--'); 
legend({'y', 'Normal range of abduction-adduction'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
ylim([-60 40]);
title('Left Hip angles: Abduction-Adduction');

figure(9)
plot(t_1000, LHipAngles_z(1:(frames/1000):frames));
hold on 
plot(t_1000, intRot_hip,'--',t_1000, extRot_hip,'--'); 
legend({'z', 'Normal range of internal-external rotation'},'Location','best');
xlabel('Time in seconds (%1000 lifting cycle)'); %tanto por mil 
ylabel('Degrees');
%ylim([-70 80]);
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

phase3_RHip_x = RHipAngles_x(idxF2:idxF3);
phase3_RHip_y = RHipAngles_y(idxF2:idxF3);
phase3_RHip_z = RHipAngles_z(idxF2:idxF3);

phase4_RHip_x = RHipAngles_x(idxF3:idxF4); 
phase4_RHip_y = RHipAngles_y(idxF3:idxF4);
phase4_RHip_z = RHipAngles_z(idxF3:idxF4);

phase5_RHip_x = RHipAngles_x(idxF4:idxF5); 
phase5_RHip_y = RHipAngles_y(idxF4:idxF5);
phase5_RHip_z = RHipAngles_z(idxF4:idxF5);

% Left Hip

phase1_LHip_x = LHipAngles_x(idxO1:idxF1);
phase1_LHip_y = LHipAngles_y(idxO1:idxF1);
phase1_LHip_z = LHipAngles_z(idxO1:idxF1);

phase2_LHip_x = LHipAngles_x(idxF1:idxF2);
phase2_LHip_y = LHipAngles_y(idxF1:idxF2);
phase2_LHip_z = LHipAngles_z(idxF1:idxF2);

phase3_LHip_x = LHipAngles_x(idxF2:idxF3);
phase3_LHip_y = LHipAngles_y(idxF2:idxF3);
phase3_LHip_z = LHipAngles_z(idxF2:idxF3);

phase4_LHip_x = LHipAngles_x(idxF3:idxF4);
phase4_LHip_y = LHipAngles_y(idxF3:idxF4);
phase4_LHip_z = LHipAngles_z(idxF3:idxF4);

phase5_LHip_x = LHipAngles_x(idxF4:idxF5);
phase5_LHip_y = LHipAngles_y(idxF4:idxF5);
phase5_LHip_z = LHipAngles_z(idxF4:idxF5);

%% time vectors definition

t2 = 0:1:(length(phase3_RHip_x)-1);
t3 = 0:1:(length(phase1_RHip_x)-1);
t4 = 0:1:(length(phase2_RHip_x)-1);
t5 = 0:1:(length(phase4_RHip_x)-1);
t6 = 0:1:(length(phase5_RHip_x)-1);

% Right Hip figures: MOVEMENT SEGMENTED INTO THE 4 PHASES

figure(10)
subplot(5,1,1)
plot (t3, phase1_RHip_x, t3, phase1_RHip_y, '--',t3,phase1_RHip_z, '.');
legend('x','y','z')
xlabel('Time in frames'); 
ylabel('Degrees');
title('RIGHT HIP: Phase 1');   % Picking up the box sagittal plane
subplot(5,1,2)
plot(t4, phase2_RHip_x,t4,phase2_RHip_y, '--',t4,phase2_RHip_z, '.');
legend('x','y','z')
xlabel('Time in frames');
ylabel('Degrees');
title('Phase 2');  % Taking steps from sagittal to frontal
subplot(5,1,3)
plot(t2, phase3_RHip_x, t2, phase3_RHip_y,'--', t2, phase3_RHip_z, '.');
legend('x','y','z')
xlabel('Time in frames'); 
ylabel('Degrees');
title('Phase 3');  % Placing the box frontal plane
subplot(5,1,4)
plot(t5, phase4_RHip_x, t5, phase4_RHip_y,'--', t5, phase4_RHip_z, '.');
legend('x','y','z')
xlabel('Time in frames'); 
ylabel('Degrees');
title('Phase 4');  % Taking steps from frontal to sagittal
subplot(5,1,5)
plot(t6, phase5_RHip_x, t6, phase5_RHip_y,'--', t6, phase5_RHip_z, '.');
legend('x','y','z')
xlabel('Time in frames'); 
ylabel('Degrees');
title('Phase 5');  % Placing the box sagittal plane again

% Left Hip figures: MOVEMENT SEGMENTED INTO THE 4 PHASES

figure(11)
subplot(5,1,1)
plot (t3, phase1_LHip_x, t3, phase1_LHip_y, '--',t3,phase1_LHip_z, '.');
legend('x','y','z')
xlabel('Time in frames'); 
ylabel('Degrees');
title('LEFT HIP: Phase 1'); % Picking up the box sagittal plane
subplot(5,1,2)
plot(t4, phase2_LHip_x,t4, phase2_LHip_y, '--',t4, phase2_LHip_z, '.');
legend('x','y','z')
xlabel('Time in frames'); 
ylabel('Degrees');
title('Phase 2');  % Taking steps from sagittal to frontal
subplot(5,1,3)
plot(t2, phase3_LHip_x, t2, phase3_LHip_y,'--', t2, phase3_LHip_z, '.');
legend('x','y','z')
xlabel('Time in frames'); 
ylabel('Degrees');
title('Phase 3');  % Placing the box frontal plane
subplot(5,1,4)
plot(t5, phase4_LHip_x, t5, phase4_LHip_y,'--', t5, phase4_LHip_z, '.');
legend('x','y','z')
xlabel('Time in frames'); 
ylabel('Degrees');
title('Phase 4');  % Taking steps from frontal to sagittal
subplot(5,1,5)
plot(t6, phase5_LHip_x, t6, phase5_LHip_y,'--', t6, phase5_LHip_z, '.');
legend('x','y','z')
xlabel('Time in frames'); 
ylabel('Degrees');
title('Phase 5');  % Placing the box sagittal plane again


%% ROM calculation

% Right Hip

RHip_ph1(:,1) = phase1_RHip_x;
RHip_ph1(:,2) = phase1_RHip_y;
RHip_ph1(:,3) = phase1_RHip_z;

RHip_ph2(:,1) = phase2_RHip_x;
RHip_ph2(:,2) = phase2_RHip_y;
RHip_ph2(:,3) = phase2_RHip_z;

RHip_ph3(:,1) = phase3_RHip_x;
RHip_ph3(:,2) = phase3_RHip_y;
RHip_ph3(:,3) = phase3_RHip_z;

RHip_ph4(:,1) = phase4_RHip_x;
RHip_ph4(:,2) = phase4_RHip_y;
RHip_ph4(:,3) = phase4_RHip_z;

RHip_ph5(:,1) = phase5_RHip_x;
RHip_ph5(:,2) = phase5_RHip_y;
RHip_ph5(:,3) = phase5_RHip_z;


RHip_ph1_max = max(RHip_ph1);
RHip_ph1_min = min(RHip_ph1);

RHip_ph2_max = max(RHip_ph2);
RHip_ph2_min = min(RHip_ph2);

RHip_ph3_max = max(RHip_ph3);
RHip_ph3_min = min(RHip_ph3);

RHip_ph4_max = max(RHip_ph4);
RHip_ph4_min = min(RHip_ph4);

RHip_ph5_max = max(RHip_ph5);
RHip_ph5_min = min(RHip_ph5);

for i = 1:length(RHip_ph1_max)
    RHip_ph1_ROM(1,i) = RHip_ph1_max(i)- RHip_ph1_min(i);
end

for i = 1:length(RHip_ph1_max)
    RHip_ph2_ROM(1,i) = RHip_ph2_max(i)- RHip_ph2_min(i);
end

for i = 1:length(RHip_ph1_max)
    RHip_ph3_ROM(1,i) = RHip_ph3_max(i)- RHip_ph3_min(i);
end

for i = 1:length(RHip_ph1_max)
    RHip_ph4_ROM(1,i) = RHip_ph4_max(i)- RHip_ph4_min(i);
end

for i = 1:length(RHip_ph1_max)
    RHip_ph5_ROM(1,i) = RHip_ph5_max(i)- RHip_ph5_min(i);
end

% Left Hip

LHip_ph1(:,1) = phase1_LHip_x;
LHip_ph1(:,2) = phase1_LHip_y;
LHip_ph1(:,3) = phase1_LHip_z;

LHip_ph2(:,1) = phase2_LHip_x;
LHip_ph2(:,2) = phase2_LHip_y;
LHip_ph2(:,3) = phase2_LHip_z;

LHip_ph3(:,1) = phase3_LHip_x;
LHip_ph3(:,2) = phase3_LHip_y;
LHip_ph3(:,3) = phase3_LHip_z;

LHip_ph4(:,1) = phase4_LHip_x;
LHip_ph4(:,2) = phase4_LHip_y;
LHip_ph4(:,3) = phase4_LHip_z;

LHip_ph5(:,1) = phase5_LHip_x;
LHip_ph5(:,2) = phase5_LHip_y;
LHip_ph5(:,3) = phase5_LHip_z;


LHip_ph1_max = max(LHip_ph1);
LHip_ph1_min = min(LHip_ph1);

LHip_ph2_max = max(LHip_ph2);
LHip_ph2_min = min(LHip_ph2);

LHip_ph3_max = max(LHip_ph3);
LHip_ph3_min = min(LHip_ph3);

LHip_ph4_max = max(LHip_ph4);
LHip_ph4_min = min(LHip_ph4);

LHip_ph5_max = max(LHip_ph5);
LHip_ph5_min = min(LHip_ph5);

for i = 1:length(RHip_ph1_max)
    LHip_ph1_ROM(1,i) = LHip_ph1_max(i)- LHip_ph1_min(i);
end

for i = 1:length(RHip_ph1_max)
    LHip_ph2_ROM(1,i) = LHip_ph2_max(i)- LHip_ph2_min(i);
end

for i = 1:length(RHip_ph1_max)
    LHip_ph3_ROM(1,i) = LHip_ph3_max(i)- LHip_ph3_min(i);
end

for i = 1:length(RHip_ph1_max)
    LHip_ph4_ROM(1,i) = LHip_ph4_max(i)- LHip_ph4_min(i);
end

for i = 1:length(RHip_ph1_max)
    LHip_ph5_ROM(1,i) = LHip_ph5_max(i)- LHip_ph5_min(i);
end

% Total ROM Hip Angles

RHip_total(:,1) = RHipAngles_x;
RHip_total(:,2) = RHipAngles_y;
RHip_total(:,3) = RHipAngles_z;

RHip_total_max = max(RHip_total);
RHip_total_min = min(RHip_total);

for i = 1:length(RHip_total_max)
    RHip_total_ROM(1,i) = RHip_total_max(i)- RHip_total_min(i);
end

LHip_total(:,1) = LHipAngles_x;
LHip_total(:,2) = LHipAngles_y;
LHip_total(:,3) = LHipAngles_z;

LHip_total_max = max(LHip_total);
LHip_total_min = min(LHip_total);

for i = 1:length(LHip_total_max)
    LHip_total_ROM(1,i) = LHip_total_max(i)- LHip_total_min(i);
end

%ROM_Right_Hip = table([RHip_ph1_ROM(1);RHip_ph2_ROM(1);RHip_ph3_ROM(1);RHip_ph4_ROM(1); RHip_ph5_ROM(1)],[RHip_ph1_ROM(2);RHip_ph2_ROM(2);RHip_ph3_ROM(2);RHip_ph4_ROM(2); RHip_ph5_ROM(2)],[RHip_ph1_ROM(3);RHip_ph2_ROM(3);RHip_ph3_ROM(3);RHip_ph4_ROM(3); RHip_ph5_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Phase 1','Phase 2','Phase 3', 'Phase 4', 'Phase 5'})
ROM_Right_Hip = struct;
ROM_Right_Hip.Phase1_X=RHip_ph1_ROM(1);
ROM_Right_Hip.Phase1_Y=RHip_ph1_ROM(2);
ROM_Right_Hip.Phase1_Z=RHip_ph1_ROM(3);
ROM_Right_Hip.Phase2_X=RHip_ph2_ROM(1);
ROM_Right_Hip.Phase2_Y=RHip_ph2_ROM(2);
ROM_Right_Hip.Phase2_Z=RHip_ph2_ROM(3);
ROM_Right_Hip.Phase3_X=RHip_ph3_ROM(1);
ROM_Right_Hip.Phase3_Y=RHip_ph3_ROM(2);
ROM_Right_Hip.Phase3_Z=RHip_ph3_ROM(3);
ROM_Right_Hip.Phase4_X=RHip_ph4_ROM(1);
ROM_Right_Hip.Phase4_Y=RHip_ph4_ROM(2);
ROM_Right_Hip.Phase4_Z=RHip_ph4_ROM(3);   
ROM_Right_Hip.Phase5_X=RHip_ph5_ROM(1);
ROM_Right_Hip.Phase5_Y=RHip_ph5_ROM(2);
ROM_Right_Hip.Phase5_Z=RHip_ph5_ROM(3);

ROM_Right_Hip

%ROM_Left_Hip = table([LHip_ph1_ROM(1);LHip_ph2_ROM(1);LHip_ph3_ROM(1);LHip_ph4_ROM(1); LHip_ph5_ROM(1)],[LHip_ph1_ROM(2);LHip_ph2_ROM(2);LHip_ph3_ROM(2);LHip_ph4_ROM(2); LHip_ph5_ROM(2)],[LHip_ph1_ROM(3);LHip_ph2_ROM(3);LHip_ph3_ROM(3);LHip_ph4_ROM(3); LHip_ph5_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Phase 1','Phase 2','Phase 3', 'Phase 4', 'Phase 5'})
ROM_Left_Hip = struct;
ROM_Left_Hip.Phase1_X=LHip_ph1_ROM(1);
ROM_Left_Hip.Phase1_Y=LHip_ph1_ROM(2);
ROM_Left_Hip.Phase1_Z=LHip_ph1_ROM(3);
ROM_Left_Hip.Phase2_X=LHip_ph2_ROM(1);
ROM_Left_Hip.Phase2_Y=LHip_ph2_ROM(2);
ROM_Left_Hip.Phase2_Z=LHip_ph2_ROM(3);
ROM_Left_Hip.Phase3_X=LHip_ph3_ROM(1);
ROM_Left_Hip.Phase3_Y=LHip_ph3_ROM(2);
ROM_Left_Hip.Phase3_Z=LHip_ph3_ROM(3);
ROM_Left_Hip.Phase4_X=LHip_ph4_ROM(1);
ROM_Left_Hip.Phase4_Y=LHip_ph4_ROM(2);
ROM_Left_Hip.Phase4_Z=LHip_ph4_ROM(3);
ROM_Left_Hip.Phase5_X=LHip_ph5_ROM(1);
ROM_Left_Hip.Phase5_Y=LHip_ph5_ROM(2);
ROM_Left_Hip.Phase5_Z=LHip_ph5_ROM(3);

ROM_Left_Hip

%TotalROM_Hip = table([RHip_total_ROM(1); LHip_total_ROM(1)],[RHip_total_ROM(2); LHip_total_ROM(2)],[RHip_total_ROM(3); LHip_total_ROM(3)],'VariableNames',{'X','Y','Z'},'RowNames',{'Right','Left'})
TotalROM_Hip = struct;
TotalROM_Hip.Right_X=RHip_total_ROM(1);
TotalROM_Hip.Right_Y=RHip_total_ROM(2);
TotalROM_Hip.Right_Z=RHip_total_ROM(3);
TotalROM_Hip.Left_X=LHip_total_ROM(1);
TotalROM_Hip.Left_Y=LHip_total_ROM(2);
TotalROM_Hip.Left_Z=LHip_total_ROM(3);

TotalROM_Hip
